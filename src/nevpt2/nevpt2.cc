//
// BAGEL - Parallel electron correlation program.
// Filename: nevpt2.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <set>
#include <src/smith/prim_op.h>
#include <src/nevpt2/nevpt2.h>
#include <src/df/dfdistt.h>
#include <src/casscf/superci.h>
#include <src/casscf/qvec.h>
#include <src/parallel/resources.h>

using namespace std;
using namespace bagel;

// Reference: C. Angeli, R. Cimiraglia, and J.-P. Malrieu, J. Chem. Phys. 117, 9138 (2002).
// Notations closely follow that in the reference


NEVPT2::NEVPT2(const shared_ptr<const PTree> input, const shared_ptr<const Geometry> g, const shared_ptr<const Reference> ref) : Method(input, g, ref) {

  casscf_ = make_shared<SuperCI>(input, g, ref);
  casscf_->compute();
  ref_ = casscf_->conv_to_ref();

  cout << endl << "  === DF-NEVPT2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = idata_->get<bool>("frozen", true);
  istate_ = idata_->get<int>("istate", 0);
  ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  // if three is a aux_basis keyword, we use that basis
  abasis_ = to_lower(idata_->get<string>("aux_basis", ""));
  norm_thresh_ = idata_->get<double>("norm_thresh", 1.0e-13);

}


void NEVPT2::compute() {

  nclosed_ = ref_->nclosed() - ncore_;
  nact_ = ref_->nact();
  nvirt_ = ref_->nvirt();

  if (nclosed_+nact_ < 1) throw runtime_error("no correlated orbitals");
  if (nact_ < 1)          throw runtime_error("no active orbitals");
  if (nvirt_ < 1)         throw runtime_error("no virtuals orbitals");

  Timer timer;

  // coefficients -- will be updated later
  shared_ptr<Matrix> ccoeff = nclosed_ ? ref_->coeff()->slice_copy(ncore_, ncore_+nclosed_) : nullptr;
  shared_ptr<Matrix> acoeff =            ref_->coeff()->slice_copy(ncore_+nclosed_, ncore_+nclosed_+nact_);
  shared_ptr<Matrix> vcoeff = nvirt_   ? ref_->coeff()->slice_copy(ncore_+nclosed_+nact_, ncore_+nclosed_+nact_+nvirt_) : nullptr;

  shared_ptr<const MatView> ocoeff = ncore_+nclosed_ ? ref_->coeff()->slice(0, ncore_+nclosed_) : nullptr;

  // obtain particle RDMs
  compute_rdm();
  // compute auxiliary RDMs
  compute_asrdm();
  // compute hole RDMs
  compute_hrdm();

  timer.tick_print("RDMs, hole RDMs, others ");

  // make canonical orbitals in closed and virtual subspaces
  vector<double> veig;
  vector<double> oeig;
  shared_ptr<Matrix> coeffall;
  // fock
  shared_ptr<const Matrix> fock;
  // core fock
  shared_ptr<const Matrix> fock_c;
  // heff_p in active
  shared_ptr<const Matrix> fock_p;
  // heff_h in active (active treated as closed)
  shared_ptr<const Matrix> fock_h;
  {
    // * core Fock operator
    shared_ptr<const Matrix> hcore = make_shared<Hcore>(geom_);
    shared_ptr<const Matrix> ofockao = nclosed_+ncore_ ? make_shared<Fock<1>>(geom_, hcore, nullptr, ocoeff, /*store*/false, /*rhf*/true) : hcore;
    // * active Fock operator
    // first make a weighted coefficient
    shared_ptr<Matrix> rdm1_mat = rdm1_->copy();
    rdm1_mat->sqrt();
    rdm1_mat->delocalize();
    auto acoeffw = make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0)) * *rdm1_mat);
    auto fockao = make_shared<Fock<1>>(geom_, ofockao, nullptr, acoeffw, /*store*/false, /*rhf*/true);
    // MO Fock
    if (nclosed_) {
      Matrix omofock(*ccoeff % *fockao * *ccoeff);
      oeig.resize(nclosed_);
      omofock.diagonalize(oeig.data());
      *ccoeff *= omofock;
    } {
      Matrix vmofock(*vcoeff % *fockao * *vcoeff);
      veig.resize(nvirt_);
      vmofock.diagonalize(veig.data());
      *vcoeff *= vmofock;
    }
    coeffall = make_shared<Matrix>(acoeff->ndim(), nclosed_+nact_+nvirt_);
    if (nclosed_)
      coeffall->copy_block(0, 0             , ccoeff->ndim(), nclosed_, ccoeff);
    coeffall->copy_block(0, nclosed_      , acoeff->ndim(), nact_   , acoeff);
    coeffall->copy_block(0, nclosed_+nact_, vcoeff->ndim(), nvirt_  , vcoeff);

    fockact_   = make_shared<Matrix>(*acoeff % *fockao  * *acoeff);
    fockact_c_ = make_shared<Matrix>(*acoeff % *ofockao * *acoeff);
    fockact_->localize();
    fockact_c_->localize();

    fock   = make_shared<Matrix>(*coeffall % *fockao  * *coeffall);
    fock_c = make_shared<Matrix>(*coeffall % *ofockao * *coeffall);

    // h'eff (only 1/2 exchange in the active space)
    auto fockao_p = make_shared<Fock<1>>(geom_, ofockao, ofockao->clone(), make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0))), /*store*/false, /*rhf*/false);
    fockact_p_ = make_shared<Matrix>(*acoeff % *fockao_p * *acoeff);
    fockact_p_->localize();
    fock_p = make_shared<Matrix>(*coeffall % *fockao_p * *coeffall);

    // h''eff (treat active orbitals as closed)
    auto fockao_h = make_shared<Fock<1>>(geom_, ofockao, nullptr, acoeff, /*store*/false, /*rhf*/true);
    fockact_h_ = make_shared<Matrix>(*acoeff % *fockao_h * *acoeff);
    fockact_h_->localize();
    fock_h = make_shared<Matrix>(*coeffall % *fockao_h * *coeffall);
  }

  // set coefficient
  ccoeff_ = ccoeff; ccoeff.reset();
  acoeff_ = acoeff; acoeff.reset();
  vcoeff_ = vcoeff; vcoeff.reset();

  timer.tick_print("Fock computation");

  // implemented in nevpt2_mat.cc
  compute_ints();
  compute_kmat();

  timer.tick_print("K matrices");

  compute_abcd();

  timer.tick_print("A, B, C, and D matrices");

  // compute transformed integrals
  shared_ptr<DFDistT> fullvi;
  shared_ptr<const Matrix> fullav;
  shared_ptr<const Matrix> fullai;
  // TODO probably we want to use JKFIT for this for consistency?
  shared_ptr<const Matrix> fullaa;
  size_t memory_size;

  {
    // first compute half transformed integrals
    shared_ptr<DFHalfDist> half, halfa;
    if (abasis_.empty()) {
      if (nclosed_)
        half = geom_->df()->compute_half_transform(ccoeff_);
      halfa = geom_->df()->compute_half_transform(acoeff_);
      // used later to determine the cache size
      memory_size = geom_->df()->block(0)->size() * 2;
      mpi__->broadcast(&memory_size, 1, 0);
    } else {
      auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
      auto cgeom = make_shared<Geometry>(*geom_, info, false);
      if (nclosed_)
        half = cgeom->df()->compute_half_transform(ccoeff_);
      halfa = cgeom->df()->compute_half_transform(acoeff_);
      // used later to determine the cache size
      memory_size = cgeom->df()->block(0)->size();
      mpi__->broadcast(&memory_size, 1, 0);
    }

    // second transform for virtual index and rearrange data
    if (nclosed_) {
      // this is now (naux, nvirt_, nclosed_), distributed by nvirt_*nclosed_. Always naux*nvirt_ block is localized to one node
      shared_ptr<DFFullDist> full = half->compute_second_transform(vcoeff_)->apply_J()->swap();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc1());
      fullvi = make_shared<DFDistT>(full, dist);
      fullvi->discard_df();
      assert(fullvi->nblocks() == 1);
    }
    {
      shared_ptr<DFFullDist> full = halfa->compute_second_transform(coeffall)->apply_J();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size());
      auto fullax_all = make_shared<DFDistT>(full, dist);
      shared_ptr<const Matrix> fullax = fullax_all->replicate();

      if (nclosed_)
        fullai = fullax->slice_copy(0, nact_*nclosed_);
      fullaa = fullax->slice_copy(nact_*nclosed_, nact_*(nclosed_+nact_));
      fullav = fullax->slice_copy(nact_*(nclosed_+nact_), nact_*(nclosed_+nact_+nvirt_));
    }
  }
  const size_t naux = geom_->naux();

  cout << "    * 3-index integral transformation done" << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // make a list of static distribution
  const int myrank = mpi__->rank();
  vector<vector<tuple<int,int,int,int>>> tasks(mpi__->size());
  // distribution of closed-closed
  if (nclosed_) {
    StaticDist ijdist(nclosed_*(nclosed_+1)/2, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0, cnt = 0; i < nclosed_; ++i)
        for (int j = i; j < nclosed_; ++j, ++cnt)
          if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
            tasks[inode].push_back(make_tuple(j, i, /*mpitags*/-1,-1));
    }
  }
  // distribution of virt-virt (cheap as both involve active indices)
  {
    StaticDist ijdist(nvirt_*(nvirt_+1)/2, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0, cnt = 0; i < nvirt_; ++i)
        for (int j = i; j < nvirt_; ++j, ++cnt)
          if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
            tasks[inode].push_back(make_tuple(j+nclosed_+nact_, i+nclosed_+nact_, /*mpitags*/-1,-1));
    }
  }
  // distribution of closed (sort of cheap)
  if (nclosed_) {
    StaticDist ijdist(nclosed_, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0; i < nclosed_; ++i)
        if (i >= ijdist.start(inode) && i < ijdist.start(inode) + ijdist.size(inode))
          tasks[inode].push_back(make_tuple(i, -1, /*mpitags*/-1,-1));
    }
  }
  // distribution of virt for S_r(-1) (sort of cheap)
  {
    StaticDist ijdist(nvirt_, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0; i < nvirt_; ++i)
        if (i >= ijdist.start(inode) && i < ijdist.start(inode) + ijdist.size(inode))
          tasks[inode].push_back(make_tuple(i+nclosed_+nact_, -1, /*mpitags*/-1,-1));
    }
  }
  {
    int nmax = 0;
    for (auto& i : tasks)
      if (nmax < i.size()) nmax = i.size();
    for (auto& i : tasks) {
      const int n = i.size();
      for (int j = 0; j != nmax-n; ++j) i.push_back(make_tuple(-1,-1,-1,-1));
    }
  }
  const int nloop = tasks[0].size();

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // TODO this is identical to code in mp2.cc. Isolate

  // start communication (n fetch behind) - n is determined by memory size
  // the data is stored in a map
  map<int, shared_ptr<Matrix>> cache;
  // pair of node and set of integers
  vector<set<int>> cachetable(mpi__->size());
  vector<int> sendreqs;

  auto cache_block = [&](const int nadd, const int ndrop) {
    assert(ndrop < nadd);
    if (ndrop >= 0) {
      for (int inode = 0; inode != mpi__->size(); ++inode) {
        const int id = get<0>(tasks[inode][ndrop]);
        const int jd = get<1>(tasks[inode][ndrop]);
        // if id and jd are no longer used in the cache, delete the element
        set<int> used;
        for (int i = ndrop+1; i <= nadd; ++i) {
          used.insert(get<0>(tasks[inode][i]));
          used.insert(get<1>(tasks[inode][i]));
        }
        if (id >= 0 && id < nclosed_ && !used.count(id)) {
          if (inode == myrank) cache.erase(id);
          cachetable[inode].erase(id);
        }
        if (jd >= 0 && jd < nclosed_ && !used.count(jd)) {
          if (inode == myrank) cache.erase(jd);
          cachetable[inode].erase(jd);
        }
      }
    }
    if (nadd < nloop) {
      // issue recv requests
      auto request_one_ = [&](const int i, const int rank) {
        if (i < 0 || i >= nclosed_) return -1;
        cachetable[rank].insert(i);
        int tag = -1;
        if (cache.find(i) == cache.end() && myrank == rank) {
          const int origin = fullvi->locate(0, i*nvirt_);
          if (origin == myrank) {
            cache[i] = fullvi->get_slice(i*nvirt_, (i+1)*nvirt_).front();
          } else {
            cache[i] = make_shared<Matrix>(naux, nvirt_, true);
            tag = mpi__->request_recv(cache[i]->data(), cache[i]->size(), origin, myrank*nclosed_+i);
          }
        }
        return tag;
      };

      // issue send requests
      auto send_one_ = [&](const int i, const int dest) {
        // see if "i" is cached at dest
        if (i < 0 || i >= nclosed_ || cachetable[dest].count(i) || fullvi->locate(0, i*nvirt_) != myrank)
          return -1;
        return mpi__->request_send(fullvi->data() + (i*nvirt_-fullvi->bstart())*naux, nvirt_*naux, dest, dest*nclosed_+i);
      };

      for (int inode = 0; inode != mpi__->size(); ++inode) {
        if (inode == myrank) {
          // recieve data from other processes
          get<2>(tasks[myrank][nadd]) = request_one_(get<0>(tasks[myrank][nadd]), myrank); // receive requests
          get<3>(tasks[myrank][nadd]) = request_one_(get<1>(tasks[myrank][nadd]), myrank);
        } else {
          // send data to other processes
          const int i = send_one_(get<0>(tasks[inode][nadd]), inode); // send requests
          if (i >= 0) sendreqs.push_back(i);
          request_one_(get<0>(tasks[inode][nadd]), inode); // update cachetable
          const int j = send_one_(get<1>(tasks[inode][nadd]), inode); // send requests
          if (j >= 0) sendreqs.push_back(j);
          request_one_(get<1>(tasks[inode][nadd]), inode);
        }
      }
    }
  };
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  const int ncache = min(memory_size/(nvirt_*nvirt_), size_t(20));
  cout << "    * ncache = " << ncache << endl;
  for (int n = 0; n != min(ncache, nloop); ++n)
    cache_block(n, -1);

  // loop over tasks
  const map<string, int> sect{{"(+0)", 0}, {"(+1)", 1}, {"(-1)", 2}, {"(+2)", 3}, {"(-2)", 4}, {"(+1)'", 5}, {"(-1)'", 6}, {"(+0)'", 7}};
  array<double,8> energy;
  fill(energy.begin(), energy.end(), 0.0);

  for (int n = 0; n != nloop; ++n) {
    // take care of data. The communication should be hidden
    if (n+ncache < nloop)
      cache_block(n+ncache, n-1);

    const int i = get<0>(tasks[myrank][n]);
    const int j = get<1>(tasks[myrank][n]);

    if (i < 0 && j < 0) {
      continue;
    } else if (i < nclosed_ && j < nclosed_ && i >= 0 && j >= 0) {
      const int ti = get<2>(tasks[myrank][n]);
      const int tj = get<3>(tasks[myrank][n]);
      if (ti >= 0) mpi__->wait(ti);
      if (tj >= 0) mpi__->wait(tj);

      shared_ptr<const Matrix> iblock = cache.at(i);
      shared_ptr<const Matrix> jblock = cache.at(j);
      const Matrix mat(*iblock % *jblock);

      // active part
      shared_ptr<const Matrix> iablock = fullai->slice_copy(i*nact_, (i+1)*nact_);
      shared_ptr<const Matrix> jablock = fullai->slice_copy(j*nact_, (j+1)*nact_);
      const Matrix mat_va(*iblock % *jablock);
      const Matrix mat_av(*iablock % *jblock);
      // hole density matrix
      const Matrix mat_vaR(mat_va * *hrdm1_);
      const Matrix mat_avR(*hrdm1_ % mat_av);
      // K' matrix
      const Matrix mat_vaKp(mat_va * *kmatp_);
      const Matrix mat_avKp(*kmatp_ % mat_av);

      // S(2)ij,rs sector
      const Matrix mat_aa(*iablock % *jablock);
      Matrix mat_aaR(nact_, nact_, true);
      Matrix mat_aaK(nact_, nact_, true);
      dgemv_("N", nact_*nact_, nact_*nact_, 1.0,  hrdm2_->data(), nact_*nact_, mat_aa.data(), 1, 0.0, mat_aaR.data(), 1);
      dgemv_("N", nact_*nact_, nact_*nact_, 1.0, kmatp2_->data(), nact_*nact_, mat_aa.data(), 1, 0.0, mat_aaK.data(), 1);
      const double norm2  = (i == j ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaR.data());
      const double denom2 = (i == j ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaK.data());
      if (norm2 > norm_thresh_)
        energy[sect.at("(+2)")] += norm2 / (-denom2/norm2 + oeig[i]+oeig[j]);

      // TODO should thread
      // S(1)ij,r sector
      double en1 = 0.0;
      for (int v = 0; v != nvirt_; ++v) {
        double norm = 0.0;
        double denom = 0.0;
        for (int a = 0; a != nact_; ++a) {
          const double va = mat_va(v, a);
          const double av = mat_av(a, v);
          const double vaR = mat_vaR(v, a);
          const double avR = mat_avR(a, v);
          const double vaK = mat_vaKp(v, a);
          const double avK = mat_avKp(a, v);
          norm  += (2.0*(va*vaR + av*avR) - av*vaR - va*avR);
          denom += (2.0*(va*vaK + av*avK) - av*vaK - va*avK);
        }
        if (norm > norm_thresh_)
          en1 += norm / (-denom/norm-veig[v]+oeig[i]+oeig[j]);
      }
      if (i == j) en1 *= 0.5;
      energy[sect.at("(+1)")] += en1;

      // S(0)ij,rs sector
      double en = 0.0;
      for (int v = 0; v != nvirt_; ++v) {
        for (int u = v+1; u < nvirt_; ++u) {
          const double vu = mat(v, u);
          const double uv = mat(u, v);
          en += 2.0*(uv*uv + vu*vu - uv*vu) / (-veig[v]+oeig[i]-veig[u]+oeig[j]);
        }
        const double vv = mat(v, v);
        en += vv*vv / (-veig[v]+oeig[i]-veig[v]+oeig[j]);
      }
      if (i != j) en *= 2.0;
      energy[sect.at("(+0)")] += en;

    } else if (i >= nclosed_+nact_ && j >= nclosed_+nact_) {
      // S(-2)rs sector
      const int iv = i-nclosed_-nact_;
      const int jv = j-nclosed_-nact_;
      shared_ptr<const Matrix> iablock = fullav->slice_copy(iv*nact_, (iv+1)*nact_);
      shared_ptr<const Matrix> jablock = fullav->slice_copy(jv*nact_, (jv+1)*nact_);
      Matrix mat_aa(*iablock % *jablock);
      Matrix mat_aaR(nact_, nact_, true);
      Matrix mat_aaK(nact_, nact_, true);
      dgemv_("N", nact_*nact_, nact_*nact_, 1.0,  rdm2_->data(), nact_*nact_, mat_aa.data(), 1, 0.0, mat_aaR.data(), 1);
      dgemv_("T", nact_*nact_, nact_*nact_, 1.0, kmat2_->data(), nact_*nact_, mat_aa.data(), 1, 0.0, mat_aaK.data(), 1);
      const double norm  = (iv == jv ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaR.data());
      const double denom = (iv == jv ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaK.data());
      if (norm > norm_thresh_)
        energy[sect.at("(-2)")] += norm / (denom/norm - veig[iv] - veig[jv]);

    } else if (i >= nclosed_+nact_ && j < 0) {
      // S(-1)r sector
      shared_ptr<Matrix> ardm3_sorted = ardm3_->clone();
      shared_ptr<Matrix> ardm2_sorted = make_shared<Matrix>(nact_*nact_*nact_, nact_, true);
      SMITH::sort_indices<1,2,0,3,    0,1,1,1>(ardm2_->data(), ardm2_sorted->data(), nact_, nact_, nact_, nact_);
      SMITH::sort_indices<1,2,0,4,3,5,0,1,1,1>(ardm3_->data(), ardm3_sorted->data(), nact_, nact_, nact_, nact_, nact_, nact_);
      const int iv = i-nclosed_-nact_;
      shared_ptr<const Matrix> rblock = fullav->slice_copy(iv*nact_, (iv+1)*nact_);
      shared_ptr<const Matrix> bac = make_shared<Matrix>(*rblock % *fullaa);
      shared_ptr<Matrix> abc = make_shared<Matrix>(nact_*nact_*nact_, 1, true);
      SMITH::sort_indices<1,0,2,0,1,1,1>(bac->data(), abc->data(), nact_, nact_, nact_);
      shared_ptr<Matrix> heff = make_shared<Matrix>(nact_, 1, true);
      for (int a = 0; a != nact_; ++a)
        heff->element(a,0) = (2.0*fock_p->element(a+nclosed_, i) - fock_c->element(a+nclosed_, i));
      const double norm = abc->dot_product(*ardm3_sorted % *abc) + 2.0*heff->dot_product(*ardm2_sorted % *abc) + heff->dot_product(*rdm1_ % *heff);
      const double denom = abc->dot_product(*amat3_ % *abc) + heff->dot_product(*bmat2_ % *abc) + heff->dot_product(*cmat2_ * *abc) + heff->dot_product(*dmat1_ % *heff);
      if (norm > norm_thresh_)
        energy[sect.at("(-1)'")] += norm / (-denom/norm - veig[iv]);

    } else if (i < nclosed_ && j < 0) {
      // (g|vi) with i fixed
      shared_ptr<const Matrix> iblock = cache.at(i);
      // (g|ai) with i fixed
      shared_ptr<const Matrix> iablock = fullai->slice_copy(i*nact_, (i+1)*nact_);
      // reordered srdm
      Matrix srdm2_p(nact_*nact_, nact_*nact_);
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(srdm2_->data(), srdm2_p.data(), nact_, nact_, nact_, nact_);

      for (int r = 0; r != nvirt_; ++r) {
        shared_ptr<const Matrix> ibr = iblock->slice_copy(r, r+1);
        shared_ptr<const Matrix> rblock = fullav->slice_copy(r*nact_, (r+1)*nact_);

        // S(-1)rs sector
        for (int s = r; s != nvirt_; ++s) {
          shared_ptr<const Matrix> ibs = iblock->slice_copy(s, s+1);
          shared_ptr<const Matrix> sblock = fullav->slice_copy(s*nact_, (s+1)*nact_);
          const Matrix mat1(*ibs % *rblock); // (vi|ar) (i, r fixed)
          const Matrix mat2(*ibr % *sblock); // (vi|as) (i, s fixed)
          const Matrix mat1R(*ibs % *rblock * *rdm1_); // (vi|ar) (i, r fixed)
          const Matrix mat2R(*ibr % *sblock * *rdm1_); // (vi|as) (i, s fixed)
          const Matrix mat1K(*ibs % *rblock * *kmat_); // (vi|ar) (i, r fixed)
          const Matrix mat2K(*ibr % *sblock * *kmat_); // (vi|as) (i, s fixed)
          const double norm  = (r == s ? 1.0 : 2.0) * (mat2R.dot_product(mat2) + mat1R.dot_product(mat1) - mat2R.dot_product(mat1));
          const double denom = (r == s ? 1.0 : 2.0) * (mat2K.dot_product(mat2) + mat1K.dot_product(mat1) - mat2K.dot_product(mat1));
          if (norm > norm_thresh_)
            energy[sect.at("(-1)")] += norm / (denom/norm + oeig[i] - veig[r] - veig[s]);
        }

        // S(0)ir sector
        const Matrix mat1(*ibr % *fullaa); // (ir|ab)  as (1,nact_*nact_)
        const Matrix mat2(*rblock % *iablock); // (ra|bi) as (nact_, nact_)
        const Matrix mat1S(mat1 * *srdm2_);
        const Matrix mat1A(mat1 * *amat2_);
        const Matrix mat1Ssym(mat1S + (mat1 ^ *srdm2_));
        const Matrix mat1Asym(mat1A + (mat1 ^ *amat2_));
              Matrix mat2Sp(nact_, nact_, true);
              Matrix mat2D (nact_, nact_, true);
        dgemv_("N", nact_*nact_, nact_*nact_, 1.0, srdm2_p.data(), nact_*nact_, mat2.data(), 1, 0.0, mat2Sp.data(), 1);
        dgemv_("N", nact_*nact_, nact_*nact_, 1.0, dmat2_->data(), nact_*nact_, mat2.data(), 1, 0.0,  mat2D.data(), 1);
        const int ir = r + nclosed_ + nact_;
        const double norm = - 2.0*mat1S.dot_product(mat1) + blas::dot_product(mat1Ssym.data(), mat1Ssym.size(), mat2.data()) + mat2Sp.dot_product(mat2)
                          + 2.0*(fock->element(ir,i) - fock_c->element(ir,i))*(fock_c->element(ir,i) + fock_h->element(ir,i))
                          + 2.0*fock_c->element(ir,i)*fock_c->element(ir,i);
        const double denom = 2.0*mat1A.dot_product(mat1) - blas::dot_product(mat1Asym.data(), mat1Asym.size(), mat2.data()) + mat2D.dot_product(mat2);
        if (norm > norm_thresh_)
          energy[sect.at("(+0)'")] += norm / (-denom/norm + oeig[i] - veig[r]);
      }

      // S(1)i sector
      shared_ptr<Matrix> ardm3_sorted = srdm3_->clone();
      shared_ptr<Matrix> ardm2_sorted = make_shared<Matrix>(nact_*nact_*nact_, nact_, true);
      SMITH::sort_indices<1,2,0,3,    0,1,1,1>(srdm2_->data(), ardm2_sorted->data(), nact_, nact_, nact_, nact_);
      SMITH::sort_indices<1,2,0,4,3,5,0,1,1,1>(srdm3_->data(), ardm3_sorted->data(), nact_, nact_, nact_, nact_, nact_, nact_);
      shared_ptr<const Matrix> bac = make_shared<Matrix>(*iablock % *fullaa);
      shared_ptr<Matrix> abc = make_shared<Matrix>(nact_*nact_*nact_, 1, true);
      SMITH::sort_indices<1,0,2,0,1,1,1>(bac->data(), abc->data(), nact_, nact_, nact_);
      shared_ptr<Matrix> heff = make_shared<Matrix>(nact_, 1, true);
      for (int a = 0; a != nact_; ++a)
        heff->element(a,0) = fock_c->element(a+nclosed_, i);
      const double norm = abc->dot_product(*ardm3_sorted % *abc) + 2.0*heff->dot_product(*ardm2_sorted % *abc) + heff->dot_product(*hrdm1_ % *heff);
      const double denom = abc->dot_product(*amat3t_ % *abc) + heff->dot_product(*bmat2t_ % *abc) + heff->dot_product(*cmat2t_ * *abc) + heff->dot_product(*dmat1t_ % *heff);
      energy[sect.at("(+1)'")] += norm / (-denom/norm + oeig[i]);
    }
  }

  // just to double check that all the communition is done
  for (auto& i : sendreqs)
    mpi__->wait(i);
  // allreduce energy contributions
  mpi__->allreduce(energy.data(), energy.size());
  energy_ = accumulate(energy.begin(), energy.end(), 0.0);

  cout << "    * assembly done" << endl << endl;
  cout << "      NEVPT2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;
  for (auto& i : sect)
    cout << "          " << setw(7) << left << i.first << right << setw(15) << setprecision(10) << energy[i.second] << endl;
  cout << endl;

  energy_ += ref_->energy();
  cout << "      NEVPT2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

}

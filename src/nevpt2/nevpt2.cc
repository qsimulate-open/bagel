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
#include <src/nevpt2/nevpt2.h>
#include <src/df/dfdistt.h>
#include <src/casscf/casbfgs.h>
#include <src/parallel/resources.h>

using namespace std;
using namespace bagel;

NEVPT2::NEVPT2(const shared_ptr<const PTree> input, const shared_ptr<const Geometry> g, const shared_ptr<const Reference> ref) : Method(input, g, ref) {

  casscf_ = make_shared<CASBFGS>(input, g, ref);
  casscf_->compute();
  ref_ = casscf_->conv_to_ref();

  cout << endl << "  === DF-NEVPT2 calculation ===" << endl << endl;

  // checks for frozen core
  const bool frozen = idata_->get<bool>("frozen", true);
  ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  // if three is a aux_basis keyword, we use that basis
  abasis_ = to_lower(idata_->get<string>("aux_basis", ""));

}


void NEVPT2::compute() {

  const size_t nclosed = ref_->nclosed() - ncore_;
  const size_t nact = ref_->nact();
  const size_t nvirt = ref_->nvirt();

  if (nclosed+nact < 1) throw runtime_error("no correlated electrons");
  if (nvirt < 1)        throw runtime_error("no virtuals orbitals");

  // coefficients
  shared_ptr<Matrix> ccoeff = ref_->coeff()->slice(ncore_, ncore_+nclosed);
  shared_ptr<Matrix> acoeff = ref_->coeff()->slice(ncore_+nclosed, ncore_+nclosed+nact);
  shared_ptr<Matrix> vcoeff = ref_->coeff()->slice(ncore_+nclosed+nact, ncore_+nclosed+nact+nvirt);
  // rdm
  shared_ptr<const Matrix> rdm1 = ref_->rdm1(/*TODO hardwired state*/0)->rdm1_mat(/*nclosed*/0);
  shared_ptr<Matrix> unit = rdm1->clone(); unit->unit();
  shared_ptr<const Matrix> hrdm1 = make_shared<Matrix>(*unit*2.0 - *rdm1);

  // Hcore
  shared_ptr<const Matrix> hcore = make_shared<Hcore>(geom_);

  // make canonical orbitals in closed and virtual subspaces
  vector<double> veig(nvirt);
  vector<double> oeig(nclosed);
  shared_ptr<const Matrix> fockact;
  {
    // * core Fock operator
    shared_ptr<const Matrix> ofockao = nclosed+ncore_ ? make_shared<const Fock<1>>(geom_, hcore, nullptr, ref_->coeff()->slice(0, ncore_+nclosed), /*store*/false, /*rhf*/true) : hcore;
    // * active Fock operator
    // first make a weighted coefficient
    shared_ptr<Matrix> acoeffw = make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0)));
    shared_ptr<Matrix> rdm1mat = rdm1->copy();
    rdm1mat->sqrt();
    *acoeffw *= *rdm1mat;
    // then make a AO density matrix
    auto fockao = make_shared<Fock<1>>(geom_, ofockao, nullptr, acoeffw, /*store*/false, /*rhf*/true);
    // MO Fock
    {
      Matrix omofock(*ccoeff % *fockao * *ccoeff);
      omofock.diagonalize(oeig.data());
      *ccoeff *= omofock;
    } {
      Matrix vmofock(*vcoeff % *fockao * *vcoeff);
      vmofock.diagonalize(veig.data());
      *vcoeff *= vmofock;
    }
    fockact = make_shared<Matrix>(*acoeff % *fockao * *acoeff);
  }


  Timer timer;
  // compute transformed integrals
  shared_ptr<DFDistT> fullvi;
  shared_ptr<Matrix> fullva;
  shared_ptr<Matrix> fullai;
  size_t memory_size;

  {
    // first compute half transformed integrals
    shared_ptr<DFHalfDist> half, halfa;
    if (abasis_.empty()) {
      half = geom_->df()->compute_half_transform(ccoeff);
      halfa = geom_->df()->compute_half_transform(acoeff);
      // used later to determine the cache size
      memory_size = half->block(0)->size() * 2;
      mpi__->broadcast(&memory_size, 1, 0);
    } else {
      auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
      auto cgeom = make_shared<Geometry>(*geom_, info, false);
      half = cgeom->df()->compute_half_transform(ccoeff);
      halfa = cgeom->df()->compute_half_transform(acoeff);
      // used later to determine the cache size
      memory_size = cgeom->df()->block(0)->size();
      mpi__->broadcast(&memory_size, 1, 0);
    }

    // second transform for virtual index and rearrange data
    {
      // this is now (naux, nvirt, nclosed), distributed by nvirt*nclosed. Always naux*nvirt block is localized to one node
      shared_ptr<DFFullDist> full = half->compute_second_transform(vcoeff)->apply_J()->swap();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc1());
      fullvi = make_shared<DFDistT>(full, dist);
    }
    {
      shared_ptr<DFFullDist> full = halfa->compute_second_transform(vcoeff)->apply_J()->swap();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size());
      auto fullvat = make_shared<DFDistT>(full, dist);
      // Matrix(naux, nvirt*nact) is replicated to each node
      fullva = fullvat->replicate();
    }
    {
      shared_ptr<DFFullDist> full = halfa->compute_second_transform(ccoeff)->apply_J();
      auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size());
      auto fullait = make_shared<DFDistT>(full, dist);
      // Matrix(naux, nact*nclosed) is replicated to each node
      fullai = fullait->replicate();
    }

    fullvi->discard_df();
  }
  assert(fullvi->nblocks() == 1);
  const size_t naux = fullvi->naux();

  cout << "    * 3-index integral transformation done" << endl;

  // make a list of static distribution of ij
  const int myrank = mpi__->rank();
  vector<vector<tuple<int,int,int,int>>> tasks(mpi__->size());
  {
    int nmax = 0;
    StaticDist ijdist(nclosed*(nclosed+1)/2, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0, cnt = 0; i < nclosed; ++i)
        for (int j = i; j < nclosed; ++j, ++cnt)
          if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
            tasks[inode].push_back(make_tuple(j, i, /*mpitags*/-1,-1));
      if (tasks[inode].size() > nmax) nmax = tasks[inode].size();
    }
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
        if (!used.count(id)) {
          if (inode == myrank) cache.erase(id);
          cachetable[inode].erase(id);
        }
        if (!used.count(jd)) {
          if (inode == myrank) cache.erase(jd);
          cachetable[inode].erase(jd);
        }
      }
    }
    if (nadd < nloop) {
      // issue recv requests
      auto request_one_ = [&](const int i, const int rank) {
        if (i < 0) return -1;
        cachetable[rank].insert(i);
        int tag = -1;
        if (cache.find(i) == cache.end() && myrank == rank) {
          const int origin = fullvi->locate(0, i*nvirt);
          if (origin == myrank) {
            cache[i] = fullvi->get_slice(i*nvirt, (i+1)*nvirt).front();
          } else {
            cache[i] = make_shared<Matrix>(naux, nvirt, true);
            tag = mpi__->request_recv(cache[i]->data(), cache[i]->size(), origin, myrank*nclosed+i);
          }
        }
        return tag;
      };

      // issue send requests
      auto send_one_ = [&](const int i, const int dest) {
        // see if "i" is cached at dest
        if (i < 0 || cachetable[dest].count(i) || fullvi->locate(0, i*nvirt) != myrank)
          return -1;
        return mpi__->request_send(fullvi->data() + (i*nvirt-fullvi->bstart())*naux, nvirt*naux, dest, dest*nclosed+i);
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

  const int ncache = min(memory_size/(nvirt*nvirt), size_t(20));
  cout << "    * ncache = " << ncache << endl;
  for (int n = 0; n != min(ncache, nloop); ++n)
    cache_block(n, -1);

  // loop over tasks
  energy_ = 0;
  double __debug__norm = 0.0;
  for (int n = 0; n != nloop; ++n) {
    // take care of data. The communication should be hidden
    if (n+ncache < nloop)
      cache_block(n+ncache, n-1);

    const int i = get<0>(tasks[myrank][n]);
    const int j = get<1>(tasks[myrank][n]);
    if (i < 0 || j < 0) continue;

    const int ti = get<2>(tasks[myrank][n]);
    const int tj = get<3>(tasks[myrank][n]);
    if (ti >= 0) mpi__->wait(ti);
    if (tj >= 0) mpi__->wait(tj);

    shared_ptr<const Matrix> iblock = cache.at(i);
    shared_ptr<const Matrix> jblock = cache.at(j);
    const Matrix mat(*iblock % *jblock);

    // active part
    shared_ptr<const Matrix> iablock = fullai->slice(i*nact, (i+1)*nact);
    shared_ptr<const Matrix> jablock = fullai->slice(j*nact, (j+1)*nact);
    const Matrix mat_va(*iblock % *jablock);
    const Matrix mat_av(*iablock % *jblock);
    // hole density matrix
    const Matrix mat_vaR(mat_va * *hrdm1);
    const Matrix mat_avR(*hrdm1 % mat_av);

    double en1 = 0.0;
    for (int v = 0; v != nvirt; ++v) {
      double norm = 0.0;
      for (int a = 0; a != nact; ++a) {
        const double va = mat_va(v, a);
        const double av = mat_av(a, v);
        const double vaR = mat_vaR(v, a);
        const double avR = mat_avR(a, v);
        norm += (2.0*(va*vaR + av*avR) - av*vaR + va*avR);
      }
      en1 += norm;
    }
    if (i != j) en1 *= 2.0;
    __debug__norm += en1;


    // TODO should thread
    double en = 0.0;
    for (int v = 0; v != nvirt; ++v) {
      for (int u = v+1; u < nvirt; ++u) {
        const double vu = mat(v, u);
        const double uv = mat(u, v);
        en += 2.0*(uv*uv + vu*vu - uv*vu) / (-veig[v]+oeig[i]-veig[u]+oeig[j]);
      }
      const double vv = mat(v, v);
      en += vv*vv / (-veig[v]+oeig[i]-veig[v]+oeig[j]);
    }
    if (i != j) en *= 2.0;
    energy_ += en;
  }

  // just to double check that all the communition is done
  for (auto& i : sendreqs)
    mpi__->wait(i);
  // allreduce energy contributions
  mpi__->allreduce(&energy_, 1);
  mpi__->allreduce(&__debug__norm, 1);
cout << setprecision(10) << __debug__norm << endl;

  cout << "    * assembly done" << endl << endl;
  cout << "      NEVPT2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;

  energy_ += ref_->energy();
  cout << "      NEVPT2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

}

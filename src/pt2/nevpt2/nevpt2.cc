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

#include <src/pt2/nevpt2/nevpt2.h>
#include <src/pt2/mp2/mp2cache.h>
#include <src/scf/hf/fock.h>
#include <src/scf/dhf/dfock.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/multi/casscf/superci.h>
#include <src/multi/zcasscf/zcashybrid.h>
#include <src/ci/zfci/relmofile.h>
#include <src/util/prim_op.h>
#include <src/util/parallel/resources.h>

using namespace std;
using namespace bagel;

// Reference: C. Angeli, R. Cimiraglia, and J.-P. Malrieu, J. Chem. Phys. 117, 9138 (2002).
// Notations closely follow that in the reference

template<typename DataType>
NEVPT2<DataType>::NEVPT2(shared_ptr<const PTree> input, shared_ptr<const Geometry> g, shared_ptr<const Reference> ref) : Method(input, g, ref) {

  const int z2 = is_same<DataType,double>::value ? 1 : 2;

  // checks for frozen core
  const bool frozen = idata_->get<bool>("frozen", true);
  istate_ = idata_->get<int>("istate", 0);
  ncore_  = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2*z2 : 0));
  if (ncore_) cout << "    * freezing " << ncore_ << " orbital" << (ncore_^1 ? "s" : "") << endl;

  // if three is a aux_basis keyword, we use that basis
  abasis_ = to_lower(idata_->get<string>("aux_basis", ""));
  norm_thresh_ = idata_->get<double>("norm_thresh", 1.0e-13);

  // starting up
  {
    init_reference();

    nclosed_ = ref_->nclosed()*z2 - ncore_;
    nact_    = ref_->nact()*z2;
    nvirt_   = ref_->nvirt()*z2;

    if (nclosed_+nact_ < 1) throw runtime_error("no correlated orbitals");
    if (nact_ < 1)          throw runtime_error("no active orbitals");
    if (nvirt_ < 1)         throw runtime_error("no virtuals orbitals");
  }

  cout << endl << "  === DF-NEVPT2 calculation ===" << endl << endl;
}


template<>
void NEVPT2<double>::init_reference() {
  auto casscf = make_shared<SuperCI>(idata_, geom_, ref_);
  casscf->compute();
  ref_ = casscf->conv_to_ref();
}


template<>
void NEVPT2<complex<double>>::init_reference() {
  auto casscf = make_shared<ZCASHybrid>(idata_, geom_, ref_);
  casscf->compute();
  ref_ = casscf->conv_to_ref();
}


template<>
shared_ptr<Matrix> NEVPT2<double>::compute_fock(shared_ptr<const Geometry> cgeom, shared_ptr<const Matrix> hcore,
                                                const MatView coeff, const double scale_exch, const double scale_coulomb) {
  return make_shared<Fock<1>>(cgeom, hcore, nullptr, coeff, /*store_half*/false, /*rhf*/true, scale_exch, scale_coulomb);
}


template<>
shared_ptr<ZMatrix> NEVPT2<complex<double>>::compute_fock(shared_ptr<const Geometry> cgeom, shared_ptr<const ZMatrix> hcore,
                                                          const ZMatView coeff, const double scale_exch, const double scale_coulomb) {
  auto ref = dynamic_pointer_cast<const RelReference>(ref_);
  return make_shared<DFock>(cgeom, hcore, coeff, ref->gaunt(), ref->breit(), /*store_half*/false, /*robust*/ref->breit(), scale_exch, scale_coulomb);
}


template<>
shared_ptr<const Matrix> NEVPT2<double>::coeff() const {
  return ref_->coeff();
}


template<>
shared_ptr<const ZMatrix> NEVPT2<complex<double>>::coeff() const {
  return dynamic_pointer_cast<const RelReference>(ref_)->relcoeff()->block_format();
}


template<>
tuple<shared_ptr<Matrix>,VectorB> NEVPT2<double>::remove_core(shared_ptr<const Matrix> in, const VectorB& o) const {
  shared_ptr<Matrix> out = in->slice_copy(ncore_, in->mdim());
  VectorB ov(in->mdim()-ncore_);
  copy(o.begin()+ncore_, o.end(), ov.begin());
  return tie(out, ov);
}


template<>
tuple<shared_ptr<ZMatrix>,VectorB> NEVPT2<complex<double>>::remove_core(shared_ptr<const ZMatrix> in, const VectorB& o) const {
  assert(in->mdim()%2 == 0 && ncore_%2 == 0);
  const int n = (in->mdim()-ncore_)/2;
  auto out = make_shared<ZMatrix>(in->ndim(), in->mdim()-ncore_, in->localized());
  out->copy_block(0, 0, in->ndim(), n, in->slice(ncore_/2, ncore_/2 + n));
  out->copy_block(0, n, in->ndim(), n, in->slice(ncore_+n, ncore_ + 2*n));

  VectorB ov(n*2);
  for (int i = 0; i != n; ++i) {
    ov(i) = o(i+ncore_/2);
    ov(i+n) = o(i+ncore_+n);
  }
  return tie(out, ov);
}


template<>
tuple<shared_ptr<DFDistT>, shared_ptr<DFDistT>>
  NEVPT2<double>::compute_full_nevpt2(shared_ptr<const Geometry> cgeom, shared_ptr<const Matrix> ccoeff, shared_ptr<const Matrix> acoeff,
                                       shared_ptr<const Matrix> vcoeff, shared_ptr<const Matrix> coeffall) const {
  shared_ptr<DFDistT> fullvi, fullax;
  if (nclosed_) {
    // this is now (naux, nvirt_, nclosed_), distributed by nvirt_*nclosed_. Always naux*nvirt_ block is localized to one node
    shared_ptr<const DFHalfDist> half = cgeom->df()->compute_half_transform(ccoeff);
    shared_ptr<DFFullDist> full = half->compute_second_transform(vcoeff)->apply_J()->swap();
    auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc1());
    fullvi = make_shared<DFDistT>(full, dist);
    fullvi->discard_df();
    assert(fullvi->nblocks() == 1);
  }
  {
    shared_ptr<const DFHalfDist> halfa = cgeom->df()->compute_half_transform(acoeff);
    shared_ptr<DFFullDist> full = halfa->compute_second_transform(coeffall)->apply_J();
    auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size());
    fullax = make_shared<DFDistT>(full, dist);
    fullax->discard_df();
    assert(fullax->nblocks() == 1);
  }
  return tie(fullvi, fullax);
}


template<>
tuple<shared_ptr<RelDFFullT>, shared_ptr<RelDFFullT>>
  NEVPT2<complex<double>>::compute_full_nevpt2(shared_ptr<const Geometry> cgeom, shared_ptr<const ZMatrix> ccoeff, shared_ptr<const ZMatrix> acoeff,
                                                shared_ptr<const ZMatrix> vcoeff, shared_ptr<const ZMatrix> coeffall) const {
  shared_ptr<RelDFFullT> fullvi, fullax;
  if (nclosed_) {
    // this is now (naux, nvirt_, nclosed_), distributed by nvirt_*nclosed_. Always naux*nvirt_ block is localized to one node
    list<shared_ptr<RelDFHalf>> half, half2;
    tie(half, half2) = RelMOFile::compute_half(cgeom, ccoeff, /*gaunt*/false, /*breit*/false);
    shared_ptr<RelDFFull> full = RelMOFile::compute_full(vcoeff, half, /*appj*/true)->swap();
    auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size(), full->nocc1());
    fullvi = make_shared<RelDFFullT>(full, dist);
    fullvi->discard_df();
  }
  {
    list<shared_ptr<RelDFHalf>> half, half2;
    tie(half, half2) = RelMOFile::compute_half(cgeom, acoeff, /*gaunt*/false, /*breit*/false);
    shared_ptr<RelDFFull> full = RelMOFile::compute_full(coeffall, half, /*appj*/true);
    auto dist = make_shared<StaticDist>(full->nocc1()*full->nocc2(), mpi__->size());
    fullax = make_shared<RelDFFullT>(full, dist);
    fullax->discard_df();
  }
  return tie(fullvi, fullax);
}


template<typename DataType>
void NEVPT2<DataType>::compute() {

  Timer timer;

  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;
  shared_ptr<const Geometry> cgeom;
  // first compute half transformed integrals
  if (abasis_.empty()) {
    cgeom = geom_;
  } else {
    auto info = make_shared<PTree>(); info->put("df_basis", abasis_);
    cgeom = make_shared<Geometry>(*geom_, info, false);
  }
  if (!is_same<DataType,double>::value) {
    auto ref = dynamic_pointer_cast<const RelReference>(ref_);
    cgeom = cgeom->relativistic(ref->gaunt());
  }

  // coefficients -- will be updated later
  shared_ptr<MatType> ccoeff = nclosed_ ? coeff()->slice_copy(0, ncore_+nclosed_) : nullptr;
  shared_ptr<MatType> vcoeff = nvirt_   ? coeff()->slice_copy(ncore_+nclosed_+nact_, ncore_+nclosed_+nact_+nvirt_) : nullptr;

  shared_ptr<const MatType> acoeff = coeff()->slice_copy(ncore_+nclosed_, ncore_+nclosed_+nact_);

  // obtain particle RDMs
  compute_rdm();
  // compute auxiliary RDMs
  compute_asrdm();
  // compute hole RDMs
  compute_hrdm();

  timer.tick_print("RDMs, hole RDMs, others ");

  /////////////////////////////////////////////////////////////////////////////////////
  // make canonical orbitals in closed and virtual subspaces
  VectorB veig(nvirt_);
  VectorB oeig(ncore_+nclosed_);
  shared_ptr<MatType> coeffall;
  // fock
  shared_ptr<const MatType> fock;
  // core fock
  shared_ptr<const MatType> fock_c;
  // heff_p in active
  shared_ptr<const MatType> fock_p;
  // heff_h in active (active treated as closed)
  shared_ptr<const MatType> fock_h;
  {
    // * core Fock operator
    shared_ptr<const MatType> hcore = make_shared<typename conditional<is_same<DataType,double>::value,Hcore,RelHcore>::type>(cgeom);
    shared_ptr<const MatType> ofockao = nclosed_+ncore_ ?  compute_fock(cgeom, hcore, coeff()->slice(0, ncore_+nclosed_)) : hcore;
    // * active Fock operator
    // first make a weighted coefficient
    MatType rdm1_mat(*rdm1_->get_conjg());
    rdm1_mat.sqrt();
    rdm1_mat.delocalize();
    const MatType acoeffw(*acoeff * (1.0/sqrt(fac2)) * rdm1_mat);
    auto fockao = compute_fock(cgeom, ofockao, acoeffw);
    // MO Fock
    if (nclosed_) {
      DiagType omofock(*ccoeff % *fockao * *ccoeff);
      omofock.diagonalize(oeig);
      *ccoeff *= omofock;
      // invoking frozen core
      if (ncore_)
        tie(ccoeff, oeig) = remove_core(ccoeff, oeig);
    } {
      DiagType vmofock(*vcoeff % *fockao * *vcoeff);
      vmofock.diagonalize(veig);
      *vcoeff *= vmofock;
    }
    coeffall = make_shared<MatType>(acoeff->ndim(), nclosed_+nact_+nvirt_);
    if (nclosed_)
      coeffall->copy_block(0, 0           , ccoeff->ndim(), nclosed_, ccoeff);
    coeffall->copy_block(0, nclosed_      , acoeff->ndim(), nact_   , acoeff);
    coeffall->copy_block(0, nclosed_+nact_, vcoeff->ndim(), nvirt_  , vcoeff);

    fockact_   = make_shared<MatType>(*acoeff % *fockao  * *acoeff);
    fockact_c_ = make_shared<MatType>(*acoeff % *ofockao * *acoeff);
    fockact_->localize();
    fockact_c_->localize();

    fock   = make_shared<MatType>(*coeffall % *fockao  * *coeffall);
    fock_c = make_shared<MatType>(*coeffall % *ofockao * *coeffall);

    // h'eff (only exchange in the active space)
    auto fockao_p = compute_fock(cgeom, ofockao, *acoeff*(1.0/sqrt(2.0)), 1.0, 0.0); // only exchange
    fockact_p_ = make_shared<MatType>(*acoeff % *fockao_p * *acoeff);
    fockact_p_->localize();
    fock_p = make_shared<MatType>(*coeffall % *fockao_p * *coeffall);

    // h''eff (treat active orbitals as closed)
    auto fockao_h = compute_fock(cgeom, ofockao, *acoeff);
    fockact_h_ = make_shared<MatType>(*acoeff % *fockao_h * *acoeff);
    fockact_h_->localize();
    fock_h = make_shared<MatType>(*coeffall % *fockao_h * *coeffall);
  }
  timer.tick_print("Fock computation");

  /////////////////////////////////////////////////////////////////////////////////////
  // compute transformed integrals
  shared_ptr<DFType> fullvi;
  shared_ptr<const MatType> fullav;
  shared_ptr<const MatType> fullai;
  // TODO probably we want to use JKFIT for this for consistency?
  shared_ptr<const MatType> fullaa;
  size_t memory_size;
  {
    // used later to determine the cache size
    memory_size = cgeom->df()->block(0)->size(); // TODO make it robust
    mpi__->broadcast(&memory_size, 1, 0);

    shared_ptr<DFType> fullax;
    tie(fullvi, fullax) = compute_full_nevpt2(cgeom, ccoeff, acoeff, vcoeff, coeffall);

    shared_ptr<const MatType> fullaxm = fullax->replicate();

    shared_ptr<MatType> tmp;
    if (nclosed_) {
      tmp = fullaxm->slice_copy(0, nact_*nclosed_);
      blas::conj_n(tmp->data(), tmp->size());
    }
    fullai = tmp;
    fullaa = fullaxm->slice_copy(nact_*nclosed_, nact_*(nclosed_+nact_));
    fullav = fullaxm->slice_copy(nact_*(nclosed_+nact_), nact_*(nclosed_+nact_+nvirt_));
  }

  {
    {
      // set 4-index integrals within the active space
      MatType ints(nact_*nact_, nact_*nact_, true);
      btas::contract(1.0, *fullaa, {0,1}, *fullaa, {0,2}, 0.0, ints, {1,2}, false, false);
      auto tmp = ints.clone();
      sort_indices<0,2,1,3,0,1,1,1>(ints.data(), tmp->data(), nact_, nact_, nact_, nact_);
      ints2_ = tmp;
    }

    {
      // set qvec
      auto rdm2c = rdm2_->copy();
      const btas::CRange<2> range(nact_, nact_*nact_*nact_);
      auto intsv = btas::make_cview(range, ints2_->storage());
      auto rdm2v = btas::make_cview(range, rdm2_->storage());
      auto qvec = make_shared<MatType>(nact_, nact_);
      btas::contract(1.0, intsv, {0,1}, rdm2v, {2,1}, 0.0, *qvec, {0,2}, false, false);
      qvec_ = qvec;
    }
  }
  cout << "    * 3-index integral transformation done" << endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // implemented in nevpt2_mat.cc
  compute_kmat();

  timer.tick_print("K matrices");

  compute_abcd();

  timer.tick_print("A, B, C, and D matrices");

  /////////////////////////////////////////////////////////////////////////////////////
  // make a list of static distribution
  vector<vector<tuple<int,int,MP2Tag<DataType>,MP2Tag<DataType>>>> tasks(mpi__->size());
  MP2Tag<DataType> dum;
  // distribution of closed-closed
  if (nclosed_) {
    StaticDist ijdist(nclosed_*(nclosed_+1)/2, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0, cnt = 0; i < nclosed_; ++i)
        for (int j = i; j < nclosed_; ++j, ++cnt)
          if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
            tasks[inode].push_back(make_tuple(j, i, /*mpitags*/dum,dum));
    }
  }
  // distribution of virt-virt (cheap as both involve active indices)
  {
    StaticDist ijdist(nvirt_*(nvirt_+1)/2, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0, cnt = 0; i < nvirt_; ++i)
        for (int j = i; j < nvirt_; ++j, ++cnt)
          if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
            tasks[inode].push_back(make_tuple(j+nclosed_+nact_, i+nclosed_+nact_, /*mpitags*/dum,dum));
    }
  }
  // distribution of closed (sort of cheap)
  if (nclosed_) {
    StaticDist ijdist(nclosed_, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0; i < nclosed_; ++i)
        if (i >= ijdist.start(inode) && i < ijdist.start(inode) + ijdist.size(inode))
          tasks[inode].push_back(make_tuple(i, -1, /*mpitags*/dum,dum));
    }
  }
  // distribution of virt for S_r(-1) (sort of cheap)
  {
    StaticDist ijdist(nvirt_, mpi__->size());
    for (int inode = 0; inode != mpi__->size(); ++inode) {
      for (int i = 0; i < nvirt_; ++i)
        if (i >= ijdist.start(inode) && i < ijdist.start(inode) + ijdist.size(inode))
          tasks[inode].push_back(make_tuple(i+nclosed_+nact_, -1, /*mpitags*/dum,dum));
    }
  }
  {
    int nmax = 0;
    for (auto& i : tasks)
      if (nmax < i.size()) nmax = i.size();
    for (auto& i : tasks) {
      const int n = i.size();
      for (int j = 0; j != nmax-n; ++j) i.push_back(make_tuple(-1,-1,dum,dum));
    }
  }

  MP2Cache_<DataType> cache(cgeom->naux(), nclosed_, nvirt_, fullvi, tasks);

  const int nloop = cache.nloop();
  const int ncache = max(3, min(static_cast<int>(memory_size/(nvirt_*nvirt_)), 20));
  cout << "    * ncache = " << ncache << endl;
  for (int n = 0; n != min(ncache, nloop); ++n)
    cache.block(n, -1);

  // loop over tasks
  const map<string, int> sect{{"(+0)", 0}, {"(+1)", 1}, {"(-1)", 2}, {"(+2)", 3}, {"(-2)", 4}, {"(+1)'", 5}, {"(-1)'", 6}, {"(+0)'", 7}};
  array<DataType,8> energy;
  fill(energy.begin(), energy.end(), static_cast<DataType>(0.0));

  // utility function that returns A^T B without conjugating A (identical to A%B for real).
  auto multiply_tn = [](const ViewType a, const ViewType b) {
    MatType out(a.mdim(), b.mdim(), true);
    btas::contract(1.0, a, {0,1}, b, {0,2}, 0.0, out, {1,2}, false, false);
    return move(out);
  };
  auto multiply_nt = [](const ViewType a, const ViewType b) {
    MatType out(a.ndim(), b.ndim(), true);
    btas::contract(1.0, a, {0,1}, b, {2,1}, 0.0, out, {0,2}, false, false);
    return move(out);
  };

  for (int n = 0; n != nloop; ++n) {
    // take care of data. The communication should be hidden
    if (n+ncache < nloop)
      cache.block(n+ncache, n-1);

    const int i = get<0>(cache.task(n));
    const int j = get<1>(cache.task(n));

    // wait till the data is available
    cache.data_wait(n);

    if (i < 0 && j < 0) {
      continue;
#if 1
    } else if (i < nclosed_ && j < nclosed_ && i >= 0 && j >= 0) {
      shared_ptr<const MatType> iblock = cache(i);
      shared_ptr<const MatType> jblock = cache(j);

      // S(0)ij,rs sector
      {
        const MatType mat = multiply_tn(*iblock, *jblock);
        DataType en = 0.0;
        for (int v = 0; v != nvirt_; ++v) {
          for (int u = v+1; u < nvirt_; ++u) {
            const DataType vu = mat(v, u);
            const DataType uv = mat(u, v);
            en += (fac2*(detail::conj(uv)*uv + detail::conj(vu)*vu) - 2.0*detail::real(detail::conj(uv)*vu)) / (-veig(v)+oeig(i)-veig(u)+oeig(j));
          }
          if (is_same<DataType,double>::value) {
            const DataType vv = mat(v, v);
            en += detail::conj(vv)*vv / (-veig(v)+oeig(i)-veig(v)+oeig(j));
          }
        }
        if (i != j) en *= 2.0;
        energy[sect.at("(+0)")] += en * (fac2*0.5); // 0.5 when relativistic
      }

      // active part
      const ViewType iablock = fullai->slice(i*nact_, (i+1)*nact_);
      const ViewType jablock = fullai->slice(j*nact_, (j+1)*nact_);
      // S(1)ij,r sector
      {
        const MatType mat_va = multiply_tn(*iblock, jablock);
        const MatType mat_av = multiply_tn(iablock, *jblock);
        // hole density matrix
        const MatType mat_vaR = multiply_nt(mat_va, *hrdm1_);
        const MatType mat_avR(*hrdm1_ * mat_av);
        // K' matrix
        const MatType mat_vaKp = multiply_nt(mat_va, *kmatp_);
        const MatType mat_avKp(*kmatp_ * mat_av);

        DataType en1 = 0.0;
        for (int v = 0; v != nvirt_; ++v) {
          DataType norm = 0.0;
          DataType denom = 0.0;
          for (int a = 0; a != nact_; ++a) {
            const DataType va = detail::conj(mat_va(v, a));
            const DataType av = detail::conj(mat_av(a, v));
            const DataType vaR = mat_vaR(v, a);
            const DataType avR = mat_avR(a, v);
            const DataType vaK = mat_vaKp(v, a);
            const DataType avK = mat_avKp(a, v);
            norm  += (fac2*(va*vaR + av*avR) - av*vaR - va*avR);
            denom += (fac2*(va*vaK + av*avK) - av*vaK - va*avK);
          }
          if (abs(norm) > norm_thresh_)
            en1 += norm / (-denom/norm-veig(v)+oeig(i)+oeig(j));
        }
        if (i == j) en1 *= 0.5;
        energy[sect.at("(+1)")] += en1;
      }

      // S(2)ij,rs sector
      const MatType mat_aa = multiply_tn(iablock, jablock);
      MatType mat_aaR(nact_, nact_, true);
      MatType mat_aaK(nact_, nact_, true);
      auto vmat_aaR = btas::group(mat_aaR,0,2);
      auto vmat_aaK = btas::group(mat_aaK,0,2);
      btas::contract(1.0, *hrdm2_,  {0,1}, btas::group(mat_aa,0,2), {1}, 0.0, vmat_aaR, {0});
      btas::contract(1.0, *kmatp2_, {0,1}, btas::group(mat_aa,0,2), {1}, 0.0, vmat_aaK, {0});
      const DataType norm2  = (i == j ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaR.data());
      const DataType denom2 = (i == j ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaK.data());
      if (abs(norm2) > norm_thresh_)
        energy[sect.at("(+2)")] += norm2 / (-denom2/norm2 + oeig(i)+oeig(j));


    } else if (i >= nclosed_+nact_ && j >= nclosed_+nact_) {
      // S(-2)rs sector
      const int iv = i-nclosed_-nact_;
      const int jv = j-nclosed_-nact_;
      const ViewType iablock = fullav->slice(iv*nact_, (iv+1)*nact_);
      const ViewType jablock = fullav->slice(jv*nact_, (jv+1)*nact_);
      MatType mat_aa = multiply_tn(iablock, jablock);
      MatType mat_aaR(nact_, nact_, true);
      MatType mat_aaK(nact_, nact_, true);
      auto vmat_aaR = btas::group(mat_aaR,0,2);
      auto vmat_aaK = btas::group(mat_aaK,0,2);
      btas::contract(1.0, *rdm2_,  {1,0}, btas::group(mat_aa,0,2), {1}, 0.0, vmat_aaR, {0});
      btas::contract(1.0, *kmat2_, {1,0}, btas::group(mat_aa,0,2), {1}, 0.0, vmat_aaK, {0});
      const DataType norm  = (iv == jv ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaR.data());
      const DataType denom = (iv == jv ? 0.5 : 1.0) * blas::dot_product(mat_aa.data(), mat_aa.size(), mat_aaK.data());
      if (abs(norm) > norm_thresh_)
        energy[sect.at("(-2)")] += norm / (denom/norm - veig[iv] - veig[jv]);

    } else if (i >= nclosed_+nact_ && j < 0) {
      // S(-1)r sector
      shared_ptr<MatType> ardm3_sorted = ardm3_->clone();
      shared_ptr<MatType> ardm2_sorted = make_shared<MatType>(nact_*nact_*nact_, nact_, true);
      sort_indices<1,2,0,3,    0,1,1,1>(ardm2_->data(), ardm2_sorted->data(), nact_, nact_, nact_, nact_);
      sort_indices<1,2,0,4,3,5,0,1,1,1>(ardm3_->data(), ardm3_sorted->data(), nact_, nact_, nact_, nact_, nact_, nact_);
      const int iv = i-nclosed_-nact_;
      const ViewType rblock = fullav->slice(iv*nact_, (iv+1)*nact_);
      const MatType bac(rblock % *fullaa);
      VecType abc(nact_*nact_*nact_);
      sort_indices<1,0,2,0,1,1,1>(bac.data(), abc.data(), nact_, nact_, nact_);
      VecType heff(nact_);
      for (int a = 0; a != nact_; ++a)
        heff(a) = detail::conj(2.0*fock_p->element(a+nclosed_, i) - fock_c->element(a+nclosed_, i));
      const DataType norm = abc % (*ardm3_sorted % abc) + heff % (2.0 * *ardm2_sorted % abc + *rdm1_ % heff);
      const DataType denom = abc % (*amat3_ % abc) + heff % (*bmat2_ % abc + *cmat2_ * abc + *dmat1_ % heff);
      if (abs(norm) > norm_thresh_)
        energy[sect.at("(-1)'")] += norm / (-denom/norm - veig[iv]);

    } else if (i < nclosed_ && j < 0) {
      // (g|vi) with i fixed
      shared_ptr<const MatType> iblock = cache(i);

      // (g|ai) with i fixed
      const ViewType iablock = fullai->slice(i*nact_, (i+1)*nact_);
      // reordered srdm
      MatType srdm2_p(nact_*nact_, nact_*nact_);
      sort_indices<0,2,1,3,0,1,1,1>(srdm2_->data(), srdm2_p.data(), nact_, nact_, nact_, nact_);

      for (int r = 0; r != nvirt_; ++r) {
        const ViewType ibr = iblock->slice(r, r+1);
        const ViewType rblock = fullav->slice(r*nact_, (r+1)*nact_);

        // S(-1)i,rs sector
        for (int s = r; s != nvirt_; ++s) {
          const ViewType ibs = iblock->slice(s, s+1);
          const ViewType sblock = fullav->slice(s*nact_, (s+1)*nact_);
          const MatType mat1 = multiply_tn(ibs, rblock); // (is|ar) (irs fixed)
          const MatType mat2 = multiply_tn(ibr, sblock); // (ir|as) (irs fixed)
          const MatType mat1R(mat1 * *rdm1_);
          const MatType mat2R(mat2 * *rdm1_);
          const MatType mat1K(mat1 * *kmat_);
          const MatType mat2K(mat2 * *kmat_);
          const DataType norm  = (r == s ? 1.0 : 2.0) * ((fac2*0.5)*(mat2R.dot_product(mat2) + mat1R.dot_product(mat1)) - detail::real(mat2R.dot_product(mat1)));
          const DataType denom = (r == s ? 1.0 : 2.0) * ((fac2*0.5)*(mat2K.dot_product(mat2) + mat1K.dot_product(mat1)) - detail::real(mat2K.dot_product(mat1)));
          if (abs(norm) > norm_thresh_)
            energy[sect.at("(-1)")] += norm / (denom/norm + oeig(i) - veig[r] - veig[s]);
        }

        // S(0)ir sector
        const MatType mat1 = multiply_tn(ibr, *fullaa); // (ir|ab)  as (1,nact_*nact_)
        const MatType mat2 = multiply_tn(rblock, iablock); // (ra|bi) as (nact_, nact_)
        const MatType mat1S(mat1 * *srdm2_);
        const MatType mat1A(mat1 * *amat2_);
        const MatType mat1Ssym(mat1S + (mat1 ^ *srdm2_));
        const MatType mat1Asym(mat1A + (mat1 ^ *amat2_));
              MatType mat2Sp(nact_, nact_, true);
              MatType mat2D (nact_, nact_, true);
        auto vmat2Sp = group(mat2Sp,0,2);
        auto vmat2D  = group(mat2D ,0,2);
        btas::contract(1.0, srdm2_p, {1,0}, btas::group(mat2,0,2), {1}, 0.0, vmat2Sp, {0});
        btas::contract(1.0, *dmat2_, {1,0}, btas::group(mat2,0,2), {1}, 0.0, vmat2D , {0});
        const int ir = r + nclosed_ + nact_;
        const DataType norm = detail::real(-fac2*mat1S.dot_product(mat1) + blas::dot_product(mat1Ssym.data(), mat1Ssym.size(), mat2.data()) + mat2Sp.dot_product(mat2))
                            + detail::real( fac2*detail::conj(fock->element(ir,i) - fock_c->element(ir,i))*(fock_c->element(ir,i) + fock_h->element(ir,i)))
                            + fac2*detail::conj(fock_c->element(ir,i))*fock_c->element(ir,i);
        const DataType denom = fac2*mat1A.dot_product(mat1) - blas::dot_product(mat1Asym.data(), mat1Asym.size(), mat2.data()) + mat2D.dot_product(mat2);
        if (abs(norm) > norm_thresh_)
          energy[sect.at("(+0)'")] += norm / (-denom/norm + oeig(i) - veig[r]);
      }

      // S(1)i sector
      shared_ptr<MatType> ardm3_sorted = srdm3_->clone();
      shared_ptr<MatType> ardm2_sorted = make_shared<MatType>(nact_*nact_*nact_, nact_, true);
      sort_indices<1,2,0,3,    0,1,1,1>(srdm2_->data(), ardm2_sorted->data(), nact_, nact_, nact_, nact_);
      sort_indices<1,2,0,4,3,5,0,1,1,1>(srdm3_->data(), ardm3_sorted->data(), nact_, nact_, nact_, nact_, nact_, nact_);
      const MatType bac(iablock % *fullaa);
      VecType abc(nact_*nact_*nact_);
      sort_indices<1,0,2,0,1,1,1>(bac.data(), abc.data(), nact_, nact_, nact_);
      VecType heff(nact_);
      for (int a = 0; a != nact_; ++a)
        heff(a) = fock_c->element(a+nclosed_, i);
      const DataType norm  = abc % (*ardm3_sorted % abc) + heff % (2.0 * *ardm2_sorted % abc + *hrdm1_->get_conjg() % heff);
      const DataType denom = abc % (*amat3t_ % abc)      + heff % (*bmat2t_ % abc + *cmat2t_ * abc + *dmat1t_ % heff);
      energy[sect.at("(+1)'")] += norm / (-denom/norm + oeig(i));
#endif
    }
  }

  // just to double check that all the communition is done
  cache.wait();
  // allreduce energy contributions
  mpi__->allreduce(energy.data(), energy.size());
  energy_ = detail::real(accumulate(energy.begin(), energy.end(), static_cast<DataType>(0.0)));

  cout << "    * assembly done" << endl << endl;
  cout << "      NEVPT2 correlation energy: " << fixed << setw(15) << setprecision(10) << energy_ << setw(10) << setprecision(2) << timer.tick() << endl << endl;
  for (auto& i : sect)
    cout << "          " << setw(7) << left << i.first << right << setw(15) << setprecision(10) << energy[i.second] << endl;
  cout << endl;

  energy_ += ref_->energy();
  cout << "      NEVPT2 total energy:       " << fixed << setw(15) << setprecision(10) << energy_ << endl << endl;

}

#define NEVPT2IMPL
#include <src/pt2/nevpt2/nevpt2_rdm.cpp>
#include <src/pt2/nevpt2/nevpt2_mat.cpp>
#undef NEVPT2IMPL

template class NEVPT2<double>;
template class NEVPT2<complex<double>>;

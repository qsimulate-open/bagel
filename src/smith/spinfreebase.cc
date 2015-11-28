//
// BAGEL - Parallel electron correlation program.
// Filename: spinfreebase.cc
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/tensor.h>
#include <src/smith/moint.h>
#include <src/smith/spinfreebase.h>
#include <src/smith/smith_util.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

template<typename DataType>
SpinFreeMethod<DataType>::SpinFreeMethod(shared_ptr<const SMITH_Info<DataType>> inf) : info_(inf) {
  static_assert(is_same<DataType,double>::value or is_same<DataType,complex<double>>::value,
                "illegal DataType for SpinFreeMethod");

  { // initializing madness world
    assert(!madness::initialized());
    int czero = 0;
    char** cnull;
    madness::initialize(czero, cnull);
  }

  Timer timer;
  const int max = info_->maxtile();
  if (info_->ncore() > info_->nclosed())
    throw runtime_error("frozen core has been specified but there are not enough closed orbitals");

  const int ncore2 = info_->ncore()*(is_same<DataType,double>::value ? 1 : 2);

  closed_ = IndexRange("c", info_->nclosed()-info_->ncore(), max, 0, info_->ncore());
  if (is_same<DataType,complex<double>>::value)
    closed_.merge(IndexRange("c", info_->nclosed()-info_->ncore(), max, closed_.nblock(), ncore2+closed_.size(), info_->ncore()));

  active_ = IndexRange("x", info_->nact(), min(max,10), closed_.nblock(), ncore2+closed_.size());
  if (is_same<DataType,complex<double>>::value)
    active_.merge(IndexRange("x", info_->nact(), min(max,10), closed_.nblock()+active_.nblock(), ncore2+closed_.size()+active_.size(),
                                                                                                 ncore2+closed_.size()));

  virt_ = IndexRange("a", info_->nvirt(), max, closed_.nblock()+active_.nblock(), ncore2+closed_.size()+active_.size());
  if (is_same<DataType,complex<double>>::value)
    virt_.merge(IndexRange("a", info_->nvirt(), max, closed_.nblock()+active_.nblock()+virt_.nblock(), ncore2+closed_.size()+active_.size()+virt_.size(),
                                                                                                       ncore2+closed_.size()+active_.size()));

  all_    = closed_; all_.merge(active_); all_.merge(virt_);

  // IndexRange for orbital update
  const int nact2 = info_->nact()*(is_same<DataType,double>::value ? 1 : 2);
  ortho1_ = IndexRange("o", nact2, max);
  ortho2_ = IndexRange("o", nact2*nact2, max);
  ortho3_ = IndexRange("o", nact2*nact2*nact2, max);
  ortho2t_ = IndexRange("o", nact2*nact2*(is_same<DataType,double>::value ? 2 : 1), max); // for XXCA

  rclosed_ = make_shared<const IndexRange>(closed_);
  ractive_ = make_shared<const IndexRange>(active_);
  rvirt_   = make_shared<const IndexRange>(virt_);

  // only for gradient computation
  if (info_->ciwfn() && info_->grad()) {
    // length of the ci expansion
    const size_t ci_size = info_->ref()->civectors()->data(info_->target())->size();
    ci_ = IndexRange("ci", ci_size, max);
    rci_ = make_shared<const IndexRange>(ci_);
  }

  // f1 tensor.
  {
    MOFock<DataType> fock(info_, {all_, all_});
    f1_ = fock.tensor();
    h1_ = fock.h1();
    eig_ = fock.eig();
    core_energy_ = fock.core_energy();
    // canonical orbitals within closed and virtual subspaces
    coeff_ = fock.coeff();
  }

  // v2 tensor.
  {
    IndexRange occ(closed_);  occ.merge(active_);
    IndexRange virt(active_); virt.merge(virt_);

    // in the case of MRCI, we need to include all sectors
    if (to_upper(info_->method()) == "MRCI") {
      occ = all_;
      virt = all_;
    }
    K2ext<DataType> v2k(info_, coeff_, {occ, virt, occ, virt});
    v2_ = v2k.tensor();
  }
  timer.tick_print("MO integral evaluation");

  // rdms.
  if (info_->ciwfn()) {
    auto fockact = make_shared<MatType>(active_.size(), active_.size());
    const int nclosed2 = info_->nclosed() * (is_same<DataType,double>::value ? 1 : 2);
    for (auto& i1 : active_)
      for (auto& i0 : active_) {
        auto loc = f1_->get_local({i0, i1});
        if (loc.first)
          fockact->copy_block(i0.offset()-nclosed2, i1.offset()-nclosed2, i0.size(), i1.size(), loc.second->get().begin());
      }
    fockact->allreduce();
    fockact = fockact->get_conjg();

    feed_rdm_denom(fockact);
    timer.tick_print("RDM + denominator evaluation");

    // rdm ci derivatives. Only for gradient computations
    if (info_->grad()) {
      feed_rdm_deriv(fockact);
      timer.tick_print("RDM derivative evaluation");
    }
  }

  // set e0
  e0_ = compute_e0();
}


template<typename DataType>
SpinFreeMethod<DataType>::~SpinFreeMethod() {
  madness::finalize();
}


template<>
void SpinFreeMethod<double>::feed_rdm_denom(shared_ptr<const Matrix> fockact) {
  const int nclo = info_->nclosed();
  const int nstates = info_->ciwfn()->nstates();
  rdm0all_ = make_shared<Vec<TATensor<double,0>>>();
  rdm1all_ = make_shared<Vec<TATensor<double,2>>>();
  rdm2all_ = make_shared<Vec<TATensor<double,4>>>();
  rdm3all_ = make_shared<Vec<TATensor<double,6>>>();
  rdm4all_ = make_shared<Vec<TATensor<double,8>>>();

  auto denom = make_shared<Denom<double>>(fockact, nstates, /*thresh*/1.0e-9);

  // TODO this can be reduced to half by bra-ket symmetry
  for (int ist = 0; ist != nstates; ++ist) {
    for (int jst = 0; jst != nstates; ++jst) {

      auto rdm0t = make_shared<TATensor<double,0>>(vector<IndexRange>());
      auto rdm1t = make_shared<TATensor<double,2>>(vector<IndexRange>(2,active_));
      auto rdm2t = make_shared<TATensor<double,4>>(vector<IndexRange>(4,active_));
      auto rdm3t = make_shared<TATensor<double,6>>(vector<IndexRange>(6,active_));
      auto rdm4t = make_shared<TATensor<double,8>>(vector<IndexRange>(8,active_));

      shared_ptr<const RDM<1>> rdm1;
      shared_ptr<const RDM<2>> rdm2;
      shared_ptr<const RDM<3>> rdm3;
      shared_ptr<const RDM<4>> rdm4; // TODO to be removed
      shared_ptr<const RDM<3>> frdm4;
      tie(rdm1, rdm2) = info_->rdm12(jst, ist);
      tie(rdm3, rdm4)  = info_->rdm34(jst, ist);
      tie(ignore, frdm4) = info_->rdm34f(jst, ist, fockact);

      (*rdm0t)("") = jst == ist ? 1.0 : 0.0;
      fill_block<2,double>(rdm1t, rdm1, vector<int>(2,nclo));
      fill_block<4,double>(rdm2t, rdm2, vector<int>(4,nclo));
      fill_block<6,double>(rdm3t, rdm3, vector<int>(6,nclo));
      fill_block<8,double>(rdm4t, rdm4, vector<int>(8,nclo));

      rdm0all_->emplace(jst, ist, rdm0t);
      rdm1all_->emplace(jst, ist, rdm1t);
      rdm2all_->emplace(jst, ist, rdm2t);
      rdm3all_->emplace(jst, ist, rdm3t);
      rdm4all_->emplace(jst, ist, rdm4t);

      denom->append(jst, ist, rdm1, rdm2, rdm3, frdm4);
    }
  }
  denom->compute();
  denom_ = denom;
}


template<>
void SpinFreeMethod<complex<double>>::feed_rdm_denom(shared_ptr<const ZMatrix> fockact) {
  const int nclo = info_->nclosed();
  const int nstates = info_->ciwfn()->nstates();
  rdm0all_ = make_shared<Vec<TATensor<complex<double>,0>>>();
  rdm1all_ = make_shared<Vec<TATensor<complex<double>,2>>>();
  rdm2all_ = make_shared<Vec<TATensor<complex<double>,4>>>();
  rdm3all_ = make_shared<Vec<TATensor<complex<double>,6>>>();
  rdm4all_ = make_shared<Vec<TATensor<complex<double>,8>>>();

  auto denom = make_shared<Denom<complex<double>>>(fockact, nstates, /*thresh*/1.0e-9);

  // TODO TODO not implemented proplerly yet
#if 0
  // TODO this can be reduced to half by bra-ket symmetry
  for (int ist = 0; ist != nstates; ++ist) {
    for (int jst = 0; jst != nstates; ++jst) {

      shared_ptr<const Kramers<2,ZRDM<1>>> rdm1;
      shared_ptr<const Kramers<4,ZRDM<2>>> rdm2;
      shared_ptr<const Kramers<6,ZRDM<3>>> rdm3;
      shared_ptr<const Kramers<8,ZRDM<4>>> rdm4;
      tie(rdm1, rdm2) = info_->rdm12(jst, ist);
      tie(rdm3, rdm4) = info_->rdm34(jst, ist);

      auto rdm1ex  = expand_kramers(rdm1, info_->nact());
      auto rdm2ex  = expand_kramers(rdm2, info_->nact());
      auto rdm3ex  = expand_kramers(rdm3, info_->nact());
      denom->append(jst, ist, rdm1ex, rdm2ex, rdm3ex, rdm4);

      auto rdm0t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>());
      unique_ptr<complex<double>[]> data0(new complex<double>[1]);
      data0[0] = jst == ist ? 1.0 : 0.0;
      rdm0t->put_block(data0);

      const int n = info_->nact();

#ifndef ALL_KRAMERS
      auto rdm1x = rdm1ex->clone();
      auto rdm2x = rdm2ex->clone();
      auto rdm3x = rdm3ex->clone();
      sort_indices<1,0,0,1,1,1>        (rdm1ex->data(), rdm1x->data(), 2*n, 2*n);
      sort_indices<1,0,3,2,0,1,1,1>    (rdm2ex->data(), rdm2x->data(), 2*n, 2*n, 2*n, 2*n);
      sort_indices<1,0,3,2,5,4,0,1,1,1>(rdm3ex->data(), rdm3x->data(), 2*n, 2*n, 2*n, 2*n, 2*n, 2*n);
      auto rdm1t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>(2,active_));
      auto rdm2t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>(4,active_));
      auto rdm3t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>(6,active_));
      fill_block<2,complex<double>>(rdm1t, rdm1x, vector<int>(2,nclo*2), vector<IndexRange>(2,active_));
      fill_block<4,complex<double>>(rdm2t, rdm2x, vector<int>(4,nclo*2), vector<IndexRange>(4,active_));
      fill_block<6,complex<double>>(rdm3t, rdm3x, vector<int>(6,nclo*2), vector<IndexRange>(6,active_));
#else
      shared_ptr<Kramers<2,ZRDM<1>>> rdm1x = rdm1->copy();
      shared_ptr<Kramers<4,ZRDM<2>>> rdm2x = rdm2->copy();
      shared_ptr<Kramers<6,ZRDM<3>>> rdm3x = rdm3->copy();
      auto j1 = rdm1x->begin();
      auto j2 = rdm2x->begin();
      auto j3 = rdm3x->begin();
      for (auto& i : *rdm1) sort_indices<1,0,0,1,1,1>        (i.second->data(), (*j1++).second->data(), n, n);
      for (auto& i : *rdm2) sort_indices<1,0,3,2,0,1,1,1>    (i.second->data(), (*j2++).second->data(), n, n, n, n);
      for (auto& i : *rdm3) sort_indices<1,0,3,2,5,4,0,1,1,1>(i.second->data(), (*j3++).second->data(), n, n, n, n, n, n);
      auto rdm1t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>(2,active_), true);
      auto rdm2t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>(4,active_), true);
      auto rdm3t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>(6,active_), true);
      fill_block<2,complex<double>,ZRDM<1>>(rdm1t, rdm1x, vector<int>(2,nclo*2), vector<IndexRange>(2,active_));
      fill_block<4,complex<double>,ZRDM<2>>(rdm2t, rdm2x, vector<int>(4,nclo*2), vector<IndexRange>(4,active_));
      fill_block<6,complex<double>,ZRDM<3>>(rdm3t, rdm3x, vector<int>(6,nclo*2), vector<IndexRange>(6,active_));
#endif
//#define RDM4_KRAMERS
#ifdef RDM4_KRAMERS
      auto rdm4x = make_shared<Kramers<8,ZRDM<4>>>();
      rdm4x->set_perm(rdm4->perm());
      for (auto& i : *rdm4) {
        shared_ptr<ZRDM<4>> data = i.second->clone();
        sort_indices<1,0,3,2,5,4,7,6,0,1,1,1>(i.second->data(), data->data(), n, n, n, n, n, n, n, n);
        rdm4x->emplace(i.first.perm({1,0,3,2,5,4,7,6}), data);
      }
      auto rdm4t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>(8,active_), true);
      fill_block<8,complex<double>,ZRDM<4>>(rdm4t, rdm4x, vector<int>(8,nclo*2), vector<IndexRange>(8,active_));
#else
      auto rdm4ex  = expand_kramers(rdm4, info_->nact());
      auto rdm4x = rdm4ex->clone();
      sort_indices<1,0,3,2,5,4,7,6,0,1,1,1>(rdm4ex->data(), rdm4x->data(), 2*n, 2*n, 2*n, 2*n, 2*n, 2*n, 2*n, 2*n);
      auto rdm4t = make_shared<Tensor_<complex<double>>>(vector<IndexRange>(8,active_));
      fill_block<8,complex<double>>(rdm4t, rdm4x, vector<int>(8,nclo*2), vector<IndexRange>(8,active_));
#endif

      rdm0all_->emplace(ist, jst, rdm0t);
      rdm1all_->emplace(ist, jst, rdm1t);
      rdm2all_->emplace(ist, jst, rdm2t);
      rdm3all_->emplace(ist, jst, rdm3t);
      rdm4all_->emplace(ist, jst, rdm4t);
    }
  }
  denom->compute();
  denom_ = denom;
#endif
}


template<>
void SpinFreeMethod<double>::feed_rdm_deriv(shared_ptr<const MatType> fockact) {
  shared_ptr<Dvec> rdm0d = make_shared<Dvec>(info_->ref()->civectors()->data(info_->target()), 1);
  shared_ptr<Dvec> rdm1d = info_->ref()->rdm1deriv(info_->target());
  shared_ptr<Dvec> rdm2d = info_->ref()->rdm2deriv(info_->target());
  shared_ptr<Dvec> rdm3d, rdm4d;
  // RDM4 is contracted a priori by the Fock operator
  tie(rdm3d, rdm4d) = info_->ref()->rdm34deriv(info_->target(), fockact);
  assert(rdm3d->ij() == rdm4d->ij());

  vector<IndexRange> o1 = {ci_};
  vector<IndexRange> o3 = {ci_, active_, active_};
  vector<IndexRange> o5 = {ci_, active_, active_, active_, active_};
  vector<IndexRange> o7 = {ci_, active_, active_, active_, active_, active_, active_};
  rdm0deriv_ = make_shared<TATensor<double,1>>(o1);
  rdm1deriv_ = make_shared<TATensor<double,3>>(o3);
  rdm2deriv_ = make_shared<TATensor<double,5>>(o5);
  rdm3deriv_ = make_shared<TATensor<double,7>>(o7);
  rdm4deriv_ = make_shared<TATensor<double,7>>(o7);

  const int nclo = info_->nclosed();
  vector<int> inpoff1(1,0);
  vector<int> inpoff3(2,nclo); inpoff3.push_back(0);
  vector<int> inpoff5(4,nclo); inpoff5.push_back(0);
  vector<int> inpoff7(6,nclo); inpoff7.push_back(0);

  const int nact = info_->nact();
  const btas::CRange<1> range1(rdm1d->extent(0)*rdm1d->extent(1));
  const btas::CRange<3> range3(rdm1d->extent(0)*rdm1d->extent(1), nact, nact);
  const btas::CRange<5> range5(rdm2d->extent(0)*rdm2d->extent(1), nact, nact, nact, nact);
  const btas::CRange<7> range7(rdm3d->extent(0)*rdm3d->extent(1), nact, nact, nact, nact, nact, nact);

  rdm0d->resize(range1);
  rdm1d->resize(range3);
  rdm2d->resize(range5);
  rdm3d->resize(range7);
  rdm4d->resize(range7);
  fill_block<1,double>(rdm0deriv_, rdm0d, inpoff1);
  fill_block<3,double>(rdm1deriv_, rdm1d, inpoff3);
  fill_block<5,double>(rdm2deriv_, rdm2d, inpoff5);
  fill_block<7,double>(rdm3deriv_, rdm3d, inpoff7);
  fill_block<7,double>(rdm4deriv_, rdm4d, inpoff7);
}


template<>
void SpinFreeMethod<complex<double>>::feed_rdm_deriv(shared_ptr<const MatType> fockact) {
  throw logic_error("SpinFreeMethod::feed_rdm_deriv is not implemented for relativistic cases.");
}


template<typename DataType>
void SpinFreeMethod<DataType>::set_rdm(const int ist, const int jst) {
  // ist is bra, jst is ket.
  // CAREFUL! the following is due to SMITH's convention (i.e., index are reversed)

  // TODO is this OK for complex cases?
  rdm0_ = rdm0all_->at(jst, ist);
  rdm1_ = rdm1all_->at(jst, ist);
  rdm2_ = rdm2all_->at(jst, ist);
  rdm3_ = rdm3all_->at(jst, ist);
  rdm4_ = rdm4all_->at(jst, ist);
}


template<typename DataType>
void SpinFreeMethod<DataType>::print_iteration() {
  cout << "      ---- iteration ----" << endl << endl;
}


template<typename DataType>
void SpinFreeMethod<DataType>::print_iteration(const int i, const double en, const double err, const double tim, const int ist) {
  cout << "     " << setw(4) << i;
  if (ist >= 0)
    cout << setw(4) << ist;
  cout << setw(15) << fixed << setprecision(8) << en << setw(15) << fixed << setprecision(8) << err
                                                     << setw(10) << fixed << setprecision(2) << tim << endl;
}


template<typename DataType>
void SpinFreeMethod<DataType>::print_iteration(const bool noconv) {
  cout << endl << "      -------------------" << endl;
  if (noconv) cout << "      *** Convergence not reached ***" << endl;
  cout << endl;
}


template<typename DataType>
double SpinFreeMethod<DataType>::compute_e0() {
  assert(!!f1_);
  const size_t nstates = info_->ciwfn()->nstates();
  DataType sum = 0.0;
  for (int ist = 0; ist != nstates; ++ist) {
    set_rdm(ist, ist);
    assert(!!rdm1_);
    sum += (*f1_)("x0,x1").dot((*rdm1_)("x0,x1")).get();
  }
  sum /= nstates;
  cout << "    - Zeroth order energy: " << setw(20) << setprecision(10) << sum << endl;
  return detail::real(sum);
}


// local function to compress the following
template<typename DataType>
void SpinFreeMethod<DataType>::loop_over(function<void(const Index&, const Index&, const Index&, const Index&)> func) const {
  for (auto& i3 : virt_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i2 : active_)
    for (auto& i0 : active_)
      for (auto& i3 : virt_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i0 : active_)
    for (auto& i3 : virt_)
      for (auto& i2 : closed_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : closed_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i2 : closed_)
        for (auto& i0 : closed_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i1 : virt_)
        for (auto& i0 : closed_) {
          func(i0, i1, i2, i3);
          func(i2, i1, i0, i3);
        }
  for (auto& i3 : active_)
    for (auto& i2 : active_)
      for (auto& i0 : active_)
        for (auto& i1 : virt_)
          func(i0, i1, i2, i3);
  for (auto& i3 : active_)
    for (auto& i1 : active_)
      for (auto& i0 : active_)
        for (auto& i2 : closed_)
          func(i0, i1, i2, i3);
}


template<typename DataType>
shared_ptr<TATensor<DataType,4>> SpinFreeMethod<DataType>::init_amplitude() const {
  shared_ptr<TATensor<DataType,4>> out = v2_->clone();
  auto put = [this, &out](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    auto local = out->get_local(vector<Index>{i0, i1, i2, i3});
    if (local.first) {
      auto it = local.second;
      const TiledArray::Range range = out->trange().make_tile_range(it.ordinal());
      typename TiledArray::Array<DataType,4>::value_type tile(range);
      std::fill_n(tile.begin(), tile.size(), 0.0);
      *it = tile;
    }
  };
  loop_over(put);
  return out;
}


template<typename DataType>
shared_ptr<TATensor<DataType,4>> SpinFreeMethod<DataType>::init_residual() const {
  shared_ptr<TATensor<DataType,4>> out = v2_->clone();
  auto put = [this, &out](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    auto local = out->get_local(vector<Index>{i2, i3, i0, i1});
    if (local.first) {
      auto it = local.second;
      const TiledArray::Range range = out->trange().make_tile_range(it.ordinal());
      typename TiledArray::Array<DataType,4>::value_type tile(range);
      std::fill_n(tile.begin(), tile.size(), 0.0);
      *it = tile;
    }
  };
  loop_over(put);
  return out;
}


template<typename DataType>
DataType SpinFreeMethod<DataType>::dot_product_transpose(shared_ptr<const MultiTATensor<DataType,4>> r, shared_ptr<const MultiTATensor<DataType,4>> t2) const {
  assert(r->nref() == t2->nref());
  DataType out = 0.0;
  for (int i = 0; i != r->nref(); ++i)
    out += detail::conj(r->fac(i)) * t2->fac(i)
         + dot_product_transpose(r->at(i), t2->at(i));
  return out;
}


template<typename DataType>
DataType SpinFreeMethod<DataType>::dot_product_transpose(shared_ptr<const TATensor<DataType,4>> r, shared_ptr<const TATensor<DataType,4>> t2) const {
  DataType out = 0.0;
  out += (*r)("c2,a3,c0,a1").dot((*t2)("c0,a1,c2,a3")).get();
  out += (*r)("x2,a3,x0,a1").dot((*t2)("x0,a1,x2,a3")).get();
  out += (*r)("c2,a3,x0,a1").dot((*t2)("x0,a1,c2,a3")).get();
  out += (*r)("c2,x3,c0,a1").dot((*t2)("c0,a1,c2,x3")).get();
  out += (*r)("c2,x3,c0,x1").dot((*t2)("c0,x1,c2,x3")).get();
  out += (*r)("x2,x3,c0,a1").dot((*t2)("c0,a1,x2,x3")).get();
  out += (*r)("c0,x3,x2,a1").dot((*t2)("x2,a1,c0,x3")).get();
  out += (*r)("x2,x3,x0,a1").dot((*t2)("x0,a1,x2,x3")).get();
  out += (*r)("c2,x3,x0,x1").dot((*t2)("x0,x1,c2,x3")).get();
  return out;
}


template<typename DataType>
shared_ptr<TATensor<DataType,4>> SpinFreeMethod<DataType>::diagonal(shared_ptr<TATensor<DataType,4>> r, shared_ptr<const TATensor<DataType,4>> t) const {
  assert(to_upper(info_->method()) == "CASPT2");
  if (!is_same<DataType,double>::value)
    throw logic_error("SpinFreeMethod<DataType>::diagonal is only correct for non-rel spin-adapted cases ");

  const int ncore = info_->ncore();
  const int nocc  = info_->nclosed() + info_->nact();
  const VecView eig = eig_;

  TATensor<DataType,4> i0(vector<IndexRange>{closed_, virt_, closed_, virt_}, true);
  i0("c2,a3,c0,a1") = (*t)("c0,a1,c2,a3")*8.0 - (*t)("c0,a3,c2,a1")*4.0;
  foreach_inplace(i0, [&](typename TATensor<DataType,4>::value_type& tile) {
    auto range = tile.range();
    auto lo = range.lobound();
    auto up = range.upbound();
    size_t n = 0;
    for (size_t i0 = lo[0]; i0 != up[0]; ++i0)
      for (size_t i1 = lo[1]; i1 != up[1]; ++i1)
        for (size_t i2 = lo[2]; i2 != up[2]; ++i2)
          for (size_t i3 = lo[3]; i3 != up[3]; ++i3)
            tile[n++] *= - eig(i3+ncore) + eig(i2+nocc) - eig(i1+ncore) + eig(i0+nocc);
  });
  (*r)("c2,a3,c0,a1") += i0("c2,a3,c0,a1");
  return r;
}

#define SPINFREEMETHOD_DETAIL
#include <src/smith/spinfreebase_update.cpp>
#undef SPINFREEMETHOD_DETAIL

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class SpinFreeMethod<double>;
template class SpinFreeMethod<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

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
#ifdef HAVE_MKL_H
    mkl_set_num_threads(1);
#endif
  }

  Timer timer;

  virt_   = info_->virt();
  active_ = info_->active();
  closed_ = info_->closed();
  all_    = info_->all();
  ci_     = info_->ci();
  ortho1_ = info_->ortho1();
  ortho2_ = info_->ortho2();
  ortho3_ = info_->ortho3();
  ortho2t_= info_->ortho2t();

  rclosed_ = make_shared<const IndexRange>(closed_);
  ractive_ = make_shared<const IndexRange>(active_);
  rvirt_   = make_shared<const IndexRange>(virt_);
  if (info_->ciwfn() && info_->grad())
    rci_ = make_shared<const IndexRange>(ci_);

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
    // TODO debug
    feed_rdm_ta();
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
#ifdef HAVE_MKL_H
  mkl_set_num_threads(info_->num_threads());
#endif
}


template<typename DataType>
void SpinFreeMethod<DataType>::feed_rdm_denom(shared_ptr<const MatType> fockact) {
  const int nclo = info_->nclosed();
  const int nstates = info_->ciwfn()->nstates();
  rdm0all_ = make_shared<Vec<TATensor<DataType,0>>>();
  rdm1all_ = make_shared<Vec<TATensor<DataType,2>>>();
  rdm2all_ = make_shared<Vec<TATensor<DataType,4>>>();
  rdm3all_ = make_shared<Vec<TATensor<DataType,6>>>();
  rdm4all_ = make_shared<Vec<TATensor<DataType,8>>>();

  const array<IndexRange,5> range{{active_, ortho1_, ortho2_, ortho3_, ortho2t_}};
  auto denom = make_shared<Denom<DataType>>(fockact, nstates, range, /*thresh*/1.0e-9);

  // TODO this can be reduced to half by bra-ket symmetry
  for (int ist = 0; ist != nstates; ++ist) {
    for (int jst = 0; jst != nstates; ++jst) {

      auto rdm0t = make_shared<TATensor<DataType,0>>(vector<IndexRange>());
      auto rdm1t = make_shared<TATensor<DataType,2>>(vector<IndexRange>(2,active_));
      auto rdm2t = make_shared<TATensor<DataType,4>>(vector<IndexRange>(4,active_));
      auto rdm3t = make_shared<TATensor<DataType,6>>(vector<IndexRange>(6,active_));
      auto rdm4t = make_shared<TATensor<DataType,8>>(vector<IndexRange>(8,active_));

      shared_ptr<const RDM<1,DataType>> rdm1;
      shared_ptr<const RDM<2,DataType>> rdm2;
      shared_ptr<const RDM<3,DataType>> rdm3;
      shared_ptr<const RDM<4,DataType>> rdm4; // TODO to be removed from here
      shared_ptr<const RDM<3,DataType>> frdm4;
      tie(rdm1, rdm2) = info_->rdm12(jst, ist);
      tie(rdm3, rdm4)  = info_->rdm34(jst, ist);
      tie(ignore, frdm4) = info_->rdm34f(jst, ist, fockact);

      denom->append(jst, ist, rdm1, rdm2, rdm3, frdm4);

      // in complex cases we reorder indices
      if (is_same<DataType,complex<double>>::value) {
        const int n = info_->nact()*2;
        auto rdm1x = rdm1->clone();
        auto rdm2x = rdm2->clone();
        auto rdm3x = rdm3->clone();
        auto rdm4x = rdm4->clone();
        sort_indices<1,0,0,1,1,1>            (rdm1->data(), rdm1x->data(), n, n);
        sort_indices<1,0,3,2,0,1,1,1>        (rdm2->data(), rdm2x->data(), n, n, n, n);
        sort_indices<1,0,3,2,5,4,0,1,1,1>    (rdm3->data(), rdm3x->data(), n, n, n, n, n, n);
        sort_indices<1,0,3,2,5,4,7,6,0,1,1,1>(rdm4->data(), rdm4x->data(), n, n, n, n, n, n, n, n);
        rdm1 = rdm1x;
        rdm2 = rdm2x;
        rdm3 = rdm3x;
        rdm4 = rdm4x;
      }

      (*rdm0t)("") = jst == ist ? 1.0 : 0.0;
      const int fac = is_same<DataType,double>::value ? 1.0 : 2.0;
      fill_block<2,DataType>(rdm1t, rdm1, vector<int>(2,nclo*fac));
      fill_block<4,DataType>(rdm2t, rdm2, vector<int>(4,nclo*fac));
      fill_block<6,DataType>(rdm3t, rdm3, vector<int>(6,nclo*fac));
      fill_block<8,DataType>(rdm4t, rdm4, vector<int>(8,nclo*fac));

      rdm0all_->emplace(jst, ist, rdm0t);
      rdm1all_->emplace(jst, ist, rdm1t);
      rdm2all_->emplace(jst, ist, rdm2t);
      rdm3all_->emplace(jst, ist, rdm3t);
      rdm4all_->emplace(jst, ist, rdm4t);
    }
  }
  denom->compute();
  denom_ = denom;
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
    if (local.first)
      out->init_tile(local.second);
  };
  loop_over(put);
  out->get_world().gop.fence();
  return out;
}


template<typename DataType>
shared_ptr<TATensor<DataType,4>> SpinFreeMethod<DataType>::init_residual() const {
  shared_ptr<TATensor<DataType,4>> out = v2_->clone();
  auto put = [&out](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    auto local = out->get_local(vector<Index>{i2, i3, i0, i1});
    if (local.first)
      out->init_tile(local.second);
  };
  loop_over(put);
  out->get_world().gop.fence();
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
  // TODO this should not be replicated here!
  out += (*r)("c2,a3,c0,a1").conj().dot((*t2)("c0,a1,c2,a3")).get();
  out += (*r)("x2,a3,x0,a1").conj().dot((*t2)("x0,a1,x2,a3")).get();
  out += (*r)("c2,a3,x0,a1").conj().dot((*t2)("x0,a1,c2,a3")).get();
  out += (*r)("c2,x3,c0,a1").conj().dot((*t2)("c0,a1,c2,x3")).get();
  out += (*r)("c2,x3,c0,x1").conj().dot((*t2)("c0,x1,c2,x3")).get();
  out += (*r)("x2,x3,c0,a1").conj().dot((*t2)("c0,a1,x2,x3")).get();
  out += (*r)("c0,x3,x2,a1").conj().dot((*t2)("x2,a1,c0,x3")).get();
  out += (*r)("x2,x3,x0,a1").conj().dot((*t2)("x0,a1,x2,x3")).get();
  out += (*r)("c2,x3,x0,x1").conj().dot((*t2)("x0,x1,c2,x3")).get();
  return out;
}



#define SPINFREEMETHOD_DETAIL
#include <src/smith/spinfreebase_update.cpp>
#include <src/smith/tardm.cpp>
#undef SPINFREEMETHOD_DETAIL

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class SpinFreeMethod<double>;
template class SpinFreeMethod<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

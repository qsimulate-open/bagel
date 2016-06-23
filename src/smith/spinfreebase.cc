//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spinfreebase.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <numeric>
#include <src/smith/moint.h>
#include <src/smith/spinfreebase.h>
#include <src/smith/smith_util.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

template<typename DataType>
SpinFreeMethod<DataType>::SpinFreeMethod(shared_ptr<const SMITH_Info<DataType>> inf) : info_(inf), info_orig_(info_) {
  static_assert(is_same<DataType,double>::value or is_same<DataType,complex<double>>::value,
                "illegal DataType for SpinFreeMethod");

  Timer timer;
  const int max = info_->maxtile();
  if (info_->ncore() > info_->nclosed())
    throw runtime_error("frozen core has been specified but there are not enough closed orbitals");

  const int ncore2 = info_->ncore()*(is_same<DataType,double>::value ? 1 : 2);

  closed_ = IndexRange(info_->nclosed()-info_->ncore(), max, 0, info_->ncore());
  if (is_same<DataType,complex<double>>::value)
    closed_.merge(IndexRange(info_->nclosed()-info_->ncore(), max, closed_.nblock(), ncore2+closed_.size(), info_->ncore()));

  active_ = IndexRange(info_->nact(), min(max,10), closed_.nblock(), ncore2+closed_.size());
  if (is_same<DataType,complex<double>>::value)
    active_.merge(IndexRange(info_->nact(), min(max,10), closed_.nblock()+active_.nblock(), ncore2+closed_.size()+active_.size(),
                                                                                            ncore2+closed_.size()));

  virt_ = IndexRange(info_->nvirt(), max, closed_.nblock()+active_.nblock(), ncore2+closed_.size()+active_.size());
  if (is_same<DataType,complex<double>>::value)
    virt_.merge(IndexRange(info_->nvirt(), max, closed_.nblock()+active_.nblock()+virt_.nblock(), ncore2+closed_.size()+active_.size()+virt_.size(),
                                                                                                  ncore2+closed_.size()+active_.size()));

  all_    = closed_; all_.merge(active_); all_.merge(virt_);

  rclosed_ = make_shared<const IndexRange>(closed_);
  ractive_ = make_shared<const IndexRange>(active_);
  rvirt_   = make_shared<const IndexRange>(virt_);

  // f1 tensor.
  {
    MOFock<DataType> fock(info_, {all_, all_});
    f1_ = fock.tensor();
    h1_ = fock.h1();
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

  auto fockact = make_shared<MatType>(active_.size(), active_.size());
  const int nclosed2 = info_->nclosed() * (is_same<DataType,double>::value ? 1 : 2);
  for (auto& i1 : active_)
    for (auto& i0 : active_)
      fockact->copy_block(i0.offset()-nclosed2, i1.offset()-nclosed2, i0.size(), i1.size(), f1_->get_block(i0, i1).get());
  fockact_ = fockact->get_conjg();

  // set Eref
  const int nstates = info_->ciwfn()->nstates();
  eref_ = make_shared<MatType>(nstates, nstates);
  for (int i = 0; i != nstates; ++i)
    eref_->element(i, i) = info_->ciwfn()->energy(i);

  if (nstates > 1 && info_->do_xms())
    rotate_xms();

  // rdms.
  if (info_->ciwfn()) {
    feed_rdm_denom();
    timer.tick_print("RDM + denominator evaluation");
  }

  // set e0all_
  compute_e0();

  mpi__->barrier();
}


template<>
void SpinFreeMethod<double>::rotate_xms() {
  assert(fockact_);
  const int nstates = info_->ciwfn()->nstates();
  Matrix fmn(nstates, nstates);

  for (int ist = 0; ist != nstates; ++ist) {
    for (int jst = 0; jst <= ist; ++jst) {
      // first compute 1RDM
      shared_ptr<const RDM<1>> rdm1;
      tie(rdm1, ignore) = info_->rdm12(jst, ist);
      // then assign the dot product: fmn=fij rdm1
      fmn(ist, jst) = blas::dot_product(fockact_->data(), fockact_->size(), rdm1->data());
      fmn(jst, ist) = fmn(ist, jst);
    }
  }

  // diagonalize fmn
  VectorB eig(nstates);
  fmn.diagonalize(eig);

  cout << endl;
  cout << "    * Extended multi-state CASPT2 (XMS-CASPT2)" << endl;
  cout << "      Rotation matrix:";
  for (int ist = 0; ist != nstates; ++ist) {
    cout << endl << "      ";
    for (int jst = 0; jst != nstates; ++jst)
      cout << setw(20) << setprecision(10) << fmn(ist, jst);
  }
  cout << endl << endl;

  // construct CIWfn
  shared_ptr<const CIWfn> ciwfn = info_->ciwfn();
  shared_ptr<const Dvec> dvec = ciwfn->civectors();
  shared_ptr<Dvec> new_dvec = dvec->clone();
  vector<shared_ptr<Civector<double>>> civecs = dvec->dvec();
  vector<shared_ptr<Civector<double>>> new_civecs = new_dvec->dvec();

  for (int jst =0; jst != nstates; ++jst) {
    for (int ist =0; ist != nstates; ++ist)
      new_civecs[jst]->ax_plus_y(fmn(ist,jst), civecs[ist]);
  }

  vector<double> energies(ciwfn->nstates());
  for (int i = 0; i != ciwfn->nstates(); ++i)
    energies[i] = ciwfn->energy(i);
  auto new_ciwfn = make_shared<CIWfn>(ciwfn->geom(), ciwfn->ncore(), ciwfn->nact(), ciwfn->nstates(),
                                      energies, new_dvec, ciwfn->det());

  // construct Reference
  auto new_ref = make_shared<Reference>(info_->geom(), make_shared<Coeff>(*info_->coeff()), info_->nclosed(), info_->nact(),
                                        info_->nvirt() + info_->nfrozenvirt(), info_->ref()->energy(), info_->ref()->rdm1(), info_->ref()->rdm2(),
                                        info_->ref()->rdm1_av(), info_->ref()->rdm2_av(), new_ciwfn);

  // construct SMITH_info
  info_ = make_shared<SMITH_Info<double>>(new_ref, info_);

  // update eref_
  eref_ = make_shared<Matrix>(fmn % (*eref_) *fmn);
  xmsmat_ = make_shared<Matrix>(move(fmn));
}


template<>
void SpinFreeMethod<complex<double>>::rotate_xms() {
  assert(fockact_);
  const int nstates = info_->ciwfn()->nstates();
  ZMatrix fmn(nstates, nstates);

  for (int ist = 0; ist != nstates; ++ist) {
    for (int jst = 0; jst <= ist; ++jst) {
      // first compute 1RDM
      shared_ptr<const Kramers<2,ZRDM<1>>> krdm1;
      tie(krdm1, ignore) = info_->rdm12(jst, ist);
      //shared_ptr<ZRDM<1>> rdm1 = expand_kramers(krdm1, krdm1->begin()->second->norb());
      shared_ptr<ZRDM<1>> rdm1 = expand_kramers(krdm1, active_.size());
      // then assign the dot product: fmn=fij rdm1
      fmn(ist, jst) = blas::dot_product(fockact_->data(), fockact_->size(), rdm1->data());
      fmn(jst, ist) = std::conj(fmn(ist, jst));
    }
  }

  // diagonalize fmn
  VectorB eig(nstates);
  fmn.diagonalize(eig);

  cout << endl;
  cout << "    * Extended multi-state CASPT2 (XMS-CASPT2)" << endl;
  cout << "      Rotation matrix:";
  for (int ist = 0; ist != nstates; ++ist) {
    cout << endl << "      ";
    for (int jst = 0; jst != nstates; ++jst)
      cout << setw(20) << setprecision(10) << fmn(ist, jst);
  }
  cout << endl << endl;

  // construct CIWfn
  // TODO:  Verify this chunk of code carefully
  shared_ptr<const RelCIWfn> ciwfn = info_->ciwfn();
  shared_ptr<const RelZDvec> dvec = ciwfn->civectors();
  shared_ptr<RelZDvec> new_dvec = dvec->clone();

  map<pair<int, int>, shared_ptr<Dvector<complex<double>>>> dvecs = dvec->dvecs();
  map<pair<int, int>, shared_ptr<Dvector<complex<double>>>> new_dvecs = new_dvec->dvecs();

//  map<pair<int, int>, vector<shared_ptr<Civector<complex<double>>>>> civecs;
//  map<pair<int, int>, vector<shared_ptr<Civector<complex<double>>>>> new_civecs;
//  for (auto& i: dvecs)
//    civecs.emplace(i.first, i.second->dvec());
//  for (auto& i: new_dvecs)
//    new_civecs.emplace(i.first, i.second->dvec());

//  for (int jst =0; jst != nstates; ++jst) {
//    for (int ist =0; ist != nstates; ++ist)
//      for (auto& i: new_civecs)
//        i.second[jst]->ax_plus_y(fmn(ist,jst), civecs.at(i.first)[ist]);
//  }
  for (auto& i: dvecs) {
    vector<shared_ptr<Civector<complex<double>>>> civecs = dvecs.at(i.first)->dvec();
    vector<shared_ptr<Civector<complex<double>>>> new_civecs = new_dvecs.at(i.first)->dvec();
    for (int jst =0; jst != nstates; ++jst)
      for (int ist =0; ist != nstates; ++ist)
        new_civecs[jst]->ax_plus_y(fmn(ist,jst), civecs[ist]);
  }

  vector<double> energies(ciwfn->nstates());
  for (int i = 0; i != ciwfn->nstates(); ++i)
    energies[i] = ciwfn->energy(i);
  auto new_ciwfn = make_shared<RelCIWfn>(ciwfn->geom(), ciwfn->ncore(), ciwfn->nact(), ciwfn->nstates(),
                                         energies, new_dvec, ciwfn->det());

  // construct Reference
  auto relref = dynamic_pointer_cast<const RelReference>(info_->ref());
  auto relcoeff = dynamic_pointer_cast<const RelCoeff_Striped>(info_->coeff());
  assert(relref && relcoeff);
  auto new_ref = make_shared<RelReference>(info_->geom(), make_shared<RelCoeff_Striped>(*relcoeff), relref->energy(),
                                           relref->nneg(), info_->nclosed(), info_->nact(), info_->nvirt() + info_->nfrozenvirt(), 
                                           info_->gaunt(), info_->breit(), /*kramers*/true,
                                           relref->rdm1_av(), relref->rdm2_av(), new_ciwfn);

  // construct SMITH_info
  info_ = make_shared<SMITH_Info<complex<double>>>(new_ref, info_);

  // update eref_
  eref_ = make_shared<ZMatrix>(fmn % (*eref_) *fmn);
  xmsmat_ = make_shared<ZMatrix>(move(fmn));
}


template<>
void SpinFreeMethod<double>::feed_rdm_denom() {
  const int nclo = info_->nclosed();
  const int nstates = info_->ciwfn()->nstates();
  rdm0all_ = make_shared<Vec<Tensor_<double>>>();
  rdm1all_ = make_shared<Vec<Tensor_<double>>>();
  rdm2all_ = make_shared<Vec<Tensor_<double>>>();
  rdm3all_ = make_shared<Vec<Tensor_<double>>>();
  rdm4all_ = make_shared<Vec<Tensor_<double>>>();

  assert(fockact_);
  auto denom = make_shared<Denom<double>>(fockact_, nstates, info_->thresh_overlap());

  // TODO this can be reduced to half by bra-ket symmetry
  for (int ist = 0; ist != nstates; ++ist) {
    for (int jst = 0; jst != nstates; ++jst) {

      shared_ptr<const RDM<1>> rdm1;
      shared_ptr<const RDM<2>> rdm2;
      shared_ptr<const RDM<3>> rdm3;
      shared_ptr<const RDM<4>> rdm4; // TODO to be removed
      tie(rdm1, rdm2) = info_->rdm12(jst, ist, (nstates > 1 && info_->do_xms()));
      tie(rdm3, rdm4) = info_->rdm34(jst, ist);
      shared_ptr<RDM<3>> frdm4 = rdm3->clone();
      auto rdm4v = group(group(*rdm4, 6,8), 0,6);
      auto frdm4v = group(*frdm4, 0, 6);
      contract(1.0, rdm4v, {0,1}, group(*fockact_,0,2), {1}, 0.0, frdm4v, {0});

      unique_ptr<double[]> data0(new double[1]);
      data0[0] = jst == ist ? 1.0 : 0.0;
      auto rdm0t = make_shared<Tensor_<double>>(vector<IndexRange>());
      rdm0t->allocate();
      if (rdm0t->is_local())
        rdm0t->put_block(data0);
      auto rdm1t = fill_block<2,double>(rdm1, vector<int>(2,nclo), vector<IndexRange>(2,active_));
      auto rdm2t = fill_block<4,double>(rdm2, vector<int>(4,nclo), vector<IndexRange>(4,active_));
      auto rdm3t = fill_block<6,double>(rdm3, vector<int>(6,nclo), vector<IndexRange>(6,active_));
      auto rdm4t = fill_block<8,double>(rdm4, vector<int>(8,nclo), vector<IndexRange>(8,active_));

      rdm0all_->emplace(jst, ist, rdm0t);
      rdm1all_->emplace(jst, ist, rdm1t);
      rdm2all_->emplace(jst, ist, rdm2t);
      rdm3all_->emplace(jst, ist, rdm3t);
      rdm4all_->emplace(jst, ist, rdm4t);

      if (!info_->sssr() || jst == ist)
        denom->append(jst, ist, rdm1, rdm2, rdm3, frdm4);
    }
  }
  denom->compute();
  denom_ = denom;
}


template<>
void SpinFreeMethod<complex<double>>::feed_rdm_denom() {
  const int nclo = info_->nclosed();
  const int nstates = info_->ciwfn()->nstates();
  rdm0all_ = make_shared<Vec<Tensor_<complex<double>>>>();
  rdm1all_ = make_shared<Vec<Tensor_<complex<double>>>>();
  rdm2all_ = make_shared<Vec<Tensor_<complex<double>>>>();
  rdm3all_ = make_shared<Vec<Tensor_<complex<double>>>>();
  rdm4all_ = make_shared<Vec<Tensor_<complex<double>>>>();

  assert(fockact_);
  auto denom = make_shared<Denom<complex<double>>>(fockact_, nstates, info_->thresh_overlap());

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
      rdm0t->allocate();
      if (rdm0t->is_local())
        rdm0t->put_block(data0);

      const int n = info_->nact();

//#define ALL_KRAMERS
#ifndef ALL_KRAMERS
      auto rdm1x = rdm1ex->clone();
      auto rdm2x = rdm2ex->clone();
      auto rdm3x = rdm3ex->clone();
      sort_indices<1,0,0,1,1,1>        (rdm1ex->data(), rdm1x->data(), 2*n, 2*n);
      sort_indices<1,0,3,2,0,1,1,1>    (rdm2ex->data(), rdm2x->data(), 2*n, 2*n, 2*n, 2*n);
      sort_indices<1,0,3,2,5,4,0,1,1,1>(rdm3ex->data(), rdm3x->data(), 2*n, 2*n, 2*n, 2*n, 2*n, 2*n);
      auto rdm1t = fill_block<2,complex<double>>(rdm1x, vector<int>(2,nclo*2), vector<IndexRange>(2,active_));
      auto rdm2t = fill_block<4,complex<double>>(rdm2x, vector<int>(4,nclo*2), vector<IndexRange>(4,active_));
      auto rdm3t = fill_block<6,complex<double>>(rdm3x, vector<int>(6,nclo*2), vector<IndexRange>(6,active_));
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
      auto rdm1t = fill_block<2,complex<double>,ZRDM<1>>(rdm1x, vector<int>(2,nclo*2), vector<IndexRange>(2,active_));
      auto rdm2t = fill_block<4,complex<double>,ZRDM<2>>(rdm2x, vector<int>(4,nclo*2), vector<IndexRange>(4,active_));
      auto rdm3t = fill_block<6,complex<double>,ZRDM<3>>(rdm3x, vector<int>(6,nclo*2), vector<IndexRange>(6,active_));
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
      auto rdm4t = fill_block<8,complex<double>,ZRDM<4>>(rdm4x, vector<int>(8,nclo*2), vector<IndexRange>(8,active_));
#else
      auto rdm4ex  = expand_kramers(rdm4, info_->nact());
      auto rdm4x = rdm4ex->clone();
      sort_indices<1,0,3,2,5,4,7,6,0,1,1,1>(rdm4ex->data(), rdm4x->data(), 2*n, 2*n, 2*n, 2*n, 2*n, 2*n, 2*n, 2*n);
      auto rdm4t = fill_block<8,complex<double>>(rdm4x, vector<int>(8,nclo*2), vector<IndexRange>(8,active_));
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
}


template<>
void SpinFreeMethod<double>::feed_rdm_deriv(const size_t offset, const size_t size) {
  tie(ci_,rci_, rdm0deriv_,rdm1deriv_,rdm2deriv_,rdm3deriv_,rdm4deriv_)
    = feed_rdm_deriv(info_, active_, fockact_, info_->target(), offset, size);
}


// static function so that it could be used by other classes (e.g., MSCASPT2)
template<>
tuple<IndexRange, shared_ptr<const IndexRange>, shared_ptr<Tensor_<double>>, shared_ptr<Tensor_<double>>,
                   shared_ptr<Tensor_<double>>, shared_ptr<Tensor_<double>>, shared_ptr<Tensor_<double>>>
  SpinFreeMethod<double>::feed_rdm_deriv(shared_ptr<const SMITH_Info<double>> info, const IndexRange& active,
                   shared_ptr<const Matrix> fockact, const int istate, const size_t offset, const size_t size) {

  auto ci = IndexRange(size, info->maxtile());
  auto rci = make_shared<const IndexRange>(ci);

  const int nact = info->nact();
  auto rdm0d = make_shared<VectorB>(size);
  copy_n(info->ref()->civectors()->data(istate)->data() + offset, size, rdm0d->data());

  auto rdm1d = make_shared<Matrix>(size, nact*nact);
  auto rdm2d = make_shared<Matrix>(size, nact*nact*nact*nact);
  {
    shared_ptr<Dvec> rdm1a = info->ref()->rdm1deriv(istate);
    shared_ptr<Dvec> rdm2a = info->ref()->rdm2deriv(istate);
    auto fill2 = [&offset, &size](shared_ptr<const Dvec> in, shared_ptr<Matrix> out) {
      assert(out->mdim() == in->ij());
      for (int i = 0; i != in->ij(); ++i)
        copy_n(in->data(i)->data() + offset, size, out->element_ptr(0,i));
    };
    fill2(rdm1a, rdm1d);
    fill2(rdm2a, rdm2d);
  }

  // RDM4 is contracted a priori by the Fock operator
  shared_ptr<Matrix> rdm3d, rdm4d;
  tie(rdm3d, rdm4d) = info->ref()->rdm34deriv(istate, fockact, offset, size);

  vector<IndexRange> o1 = {ci};
  vector<IndexRange> o3 = {ci, active, active};
  vector<IndexRange> o5 = {ci, active, active, active, active};
  vector<IndexRange> o7 = {ci, active, active, active, active, active, active};

  const int nclo = info->nclosed();
  const vector<int> inpoff1{0};
  const vector<int> inpoff3{0,nclo,nclo};
  const vector<int> inpoff5{0,nclo,nclo,nclo,nclo};
  const vector<int> inpoff7{0,nclo,nclo,nclo,nclo,nclo,nclo};

  const btas::CRange<1> range1(rdm1d->extent(0));
  const btas::CRange<3> range3(rdm1d->extent(0), nact, nact);
  const btas::CRange<5> range5(rdm2d->extent(0), nact, nact, nact, nact);
  const btas::CRange<7> range7(rdm3d->extent(0), nact, nact, nact, nact, nact, nact);

  static_pointer_cast<btas::Tensor1<double>>(rdm0d)->resize(range1);
  static_pointer_cast<btas::Tensor2<double>>(rdm1d)->resize(range3);
  static_pointer_cast<btas::Tensor2<double>>(rdm2d)->resize(range5);
  static_pointer_cast<btas::Tensor2<double>>(rdm3d)->resize(range7);
  static_pointer_cast<btas::Tensor2<double>>(rdm4d)->resize(range7);
  auto rdm0deriv = fill_block<1,double>(rdm0d, inpoff1, o1);
  auto rdm1deriv = fill_block<3,double>(rdm1d, inpoff3, o3);
  auto rdm2deriv = fill_block<5,double>(rdm2d, inpoff5, o5);
  auto rdm3deriv = fill_block<7,double>(rdm3d, inpoff7, o7);
  auto rdm4deriv = fill_block<7,double>(rdm4d, inpoff7, o7);

  return tie(ci, rci, rdm0deriv, rdm1deriv, rdm2deriv, rdm3deriv, rdm4deriv);
}


template<>
void SpinFreeMethod<complex<double>>::feed_rdm_deriv(const size_t offset, const size_t size) {
  throw logic_error("SpinFreeMethod::feed_rdm_deriv is not implemented for relativistic cases.");
}


template<>
tuple<IndexRange, shared_ptr<const IndexRange>, shared_ptr<Tensor_<complex<double>>>, shared_ptr<Tensor_<complex<double>>>,
          shared_ptr<Tensor_<complex<double>>>, shared_ptr<Tensor_<complex<double>>>, shared_ptr<Tensor_<complex<double>>>>
  SpinFreeMethod<complex<double>>::feed_rdm_deriv(shared_ptr<const SMITH_Info<complex<double>>> info, const IndexRange& active,
          shared_ptr<const ZMatrix> fockact, const int istate, const size_t offset, const size_t size) {
  throw logic_error("SpinFreeMethod::feed_rdm_deriv is not implemented for relativistic cases.");
  IndexRange d;
  shared_ptr<IndexRange> du;
  shared_ptr<Tensor_<complex<double>>> dum;
  return tie(d, du, dum, dum, dum, dum, dum);
}


template<typename DataType>
void SpinFreeMethod<DataType>::set_rdm(const int ist, const int jst) {
  // ist is bra, jst is ket.
  // CAREFUL! the following is due to SMITH's convention (i.e., index are reversed)
  rdm0_ = rdm0all_->at(jst, ist);
  rdm1_ = rdm1all_->at(jst, ist);
  rdm2_ = rdm2all_->at(jst, ist);
  rdm3_ = rdm3all_->at(jst, ist);
  rdm4_ = rdm4all_->at(jst, ist);

  // ensure that get_block calls are done after RDMs are set in every node
  mpi__->barrier();
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
void SpinFreeMethod<DataType>::compute_e0() {
  assert(!!f1_);
  const size_t nstates = info_->ciwfn()->nstates();
  e0all_.resize(nstates);
  for (int ist = 0; ist != nstates; ++ist) {
    DataType sum = 0.0;
    set_rdm(ist, ist);
    assert(!!rdm1_);
    for (auto& i1 : active_) {
      for (auto& i0 : active_) {
        if (f1_->is_local(i0, i1)) {
          const size_t size = i0.size() * i1.size();
          unique_ptr<DataType[]> fdata = f1_->get_block(i0, i1);
          unique_ptr<DataType[]> rdata = rdm1_->get_block(i0, i1);
          sum += blas::dot_product_noconj(fdata.get(), size, rdata.get());
        }
      }
    }
    mpi__->allreduce(&sum, 1);
    e0all_[ist] = detail::real(sum);
  }
  // printout
  cout << endl;
  for (int ist = 0; ist != nstates; ++ist)
    cout << "    * Zeroth order energy : state " << setw(2) << ist << setw(20) << setprecision(10) << e0all_[ist] << endl;
  cout << endl;
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
shared_ptr<Tensor_<DataType>> SpinFreeMethod<DataType>::init_amplitude() const {
  unordered_set<size_t> sparse;
  auto put = [&sparse](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    sparse.insert(generate_hash_key(i0, i1, i2, i3));
  };
  loop_over(put);
  return make_shared<Tensor_<DataType>>(v2_->indexrange(), /*kramers*/false, sparse, /*alloc*/true);
}


template<typename DataType>
shared_ptr<Tensor_<DataType>> SpinFreeMethod<DataType>::init_residual() const {
  unordered_set<size_t> sparse;
  auto put = [&sparse](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    sparse.insert(generate_hash_key(i2, i3, i0, i1));
  };
  loop_over(put);
  return make_shared<Tensor_<DataType>>(v2_->indexrange(), /*kramers*/false, sparse, /*alloc*/true);
}


template<typename DataType>
DataType SpinFreeMethod<DataType>::dot_product_transpose(shared_ptr<const MultiTensor_<DataType>> r, shared_ptr<const MultiTensor_<DataType>> t2) const {
  assert(r->nref() == t2->nref());
  DataType out = 0.0;
  for (int i = 0; i != r->nref(); ++i)
    out += detail::conj(r->fac(i)) * t2->fac(i)
         + dot_product_transpose(r->at(i), t2->at(i));
  return out;
}


template<typename DataType>
DataType SpinFreeMethod<DataType>::dot_product_transpose(shared_ptr<const Tensor_<DataType>> r, shared_ptr<const Tensor_<DataType>> t2) const {
  DataType out = 0.0;
  auto prod = [this, &r, &t2, &out](const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
    const size_t size = r->get_size(i2, i3, i0, i1);
    if (r->is_local(i2, i3, i0, i1) && size != 0) {
      unique_ptr<DataType[]> tmp0 = t2->get_block(i0, i1, i2, i3);
      unique_ptr<DataType[]> tmp1(new DataType[size]);
      sort_indices<2,3,0,1,0,1,1,1>(tmp0.get(), tmp1.get(), i0.size(), i1.size(), i2.size(), i3.size());

      out += blas::dot_product(r->get_block(i2, i3, i0, i1).get(), size, tmp1.get());
    }
  };
  loop_over(prod);
  mpi__->allreduce(&out, 1);
  return out;
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

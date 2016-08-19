//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_base.cc
// Copyright (C) 2016 Toru Shiozaki
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


#include <src/ci/fci/fci_base.h>

using namespace std;
using namespace bagel;

template<class CivecType, class DvecType>
void FCI_base<CivecType,DvecType>::compute_rdm12() {
  // Needs initialization here because we use daxpy.
  // For nstate_ == 1, rdm1_av_ = rdm1_->at(0).
  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(norb_);
    rdm2_av_ = make_shared<RDM<2>>(norb_);
  } else if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }
  // we need expanded lists
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, /*compressed=*/false, /*mute=*/true);
  cc_->set_det(detex);

  for (int i = 0; i != nstate_; ++i)
    compute_rdm12(i, i);

  // calculate state averaged RDMs
  if (nstate_ != 1) {
    for (int ist = 0; ist != nstate_; ++ist) {
      rdm1_av_->ax_plus_y(weight_[ist], rdm1_->at(ist));
      rdm2_av_->ax_plus_y(weight_[ist], rdm2_->at(ist));
    }
  } else {
    rdm1_av_ = rdm1_->at(0,0);
    rdm2_av_ = rdm2_->at(0,0);
  }

  cc_->set_det(det_);
}


template<class CivecType, class DvecType>
void FCI_base<CivecType,DvecType>::compute_rdm12(const int ist, const int jst) {
  if (det_->compress()) {
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }

  shared_ptr<CivecType> ccbra = cc_->data(ist);
  shared_ptr<CivecType> ccket = cc_->data(jst);

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  tie(rdm1, rdm2) = compute_rdm12_from_civec(ccbra, ccket);

  // setting to private members.
  rdm1_->emplace(ist, jst, rdm1);
  rdm2_->emplace(ist, jst, rdm2);

  cc_->set_det(det_);
}


template<class CivecType, class DvecType>
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
FCI_base<CivecType,DvecType>::compute_rdm12_av_from_dvec(shared_ptr<const DvecType> dbra, shared_ptr<const DvecType> dket, shared_ptr<const Determinants> o) const {

  if (o != nullptr) {
    dbra->set_det(o);
    dket->set_det(o);
  }

  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);

  assert(dbra->ij() == dket->ij() && dbra->det() == dket->det());

  for (int i = 0; i != dbra->ij(); ++i) {
    shared_ptr<RDM<1>> r1;
    shared_ptr<RDM<2>> r2;
    tie(r1, r2) = compute_rdm12_from_civec(dbra->data(i), dket->data(i));
    rdm1->ax_plus_y(weight_[i], r1);
    rdm2->ax_plus_y(weight_[i], r2);
  }

  if (o != nullptr) {
    dbra->set_det(det_);
    dket->set_det(det_);
  }

  return tie(rdm1, rdm2);
}


template<class CivecType, class DvecType>
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
FCI_base<CivecType,DvecType>::compute_rdm12_last_step(shared_ptr<const DvecType> dbra, shared_ptr<const DvecType> dket, shared_ptr<const CivecType> cibra) const {

  const int nri = cibra->asize()*cibra->lenb();
  const int ij  = norb_*norb_;

  // 1RDM c^dagger <I|\hat{E}|0>
  // 2RDM \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);
  {
    auto cibra_data = make_shared<VectorB>(nri);
    copy_n(cibra->data(), nri, cibra_data->data());

    auto dket_data = make_shared<Matrix>(nri, ij);
    for (int i = 0; i != ij; ++i)
      copy_n(dket->data(i)->data(), nri, dket_data->element_ptr(0, i));
    auto rdm1t = btas::group(*rdm1,0,2);
    btas::contract(1.0, *dket_data, {0,1}, *cibra_data, {0}, 0.0, rdm1t, {1}); 

    auto dbra_data = dket_data;
    if (dbra != dket) {
      dbra_data = make_shared<Matrix>(nri, ij);
      for (int i = 0; i != ij; ++i)
        copy_n(dbra->data(i)->data(), nri, dbra_data->element_ptr(0, i));
    }
    auto rdm2t = group(group(*rdm2, 2,4), 0,2);
    btas::contract(1.0, *dbra_data, {1,0}, *dket_data, {1,2}, 0.0, rdm2t, {0,2}); 
  }

  // when used with dist FCI, we have to allreduce
  if (is_same<CivecType,DistCivec>::value) {
    rdm1->allreduce();
    rdm2->allreduce();
  }

  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      copy_n(&rdm2->element(0,0,k,i), norb_*norb_, buf.get());
      blas::transpose(buf.get(), norb_, norb_, rdm2->element_ptr(0,0,k,i));
    }
  }

  // put in diagonal into 2RDM
  // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
  for (int i = 0; i != norb_; ++i)
    for (int k = 0; k != norb_; ++k)
      for (int j = 0; j != norb_; ++j)
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  return tie(rdm1, rdm2);
}


// note that this does not transform internal integrals (since it is not needed in CASSCF).
template<class CivecType, class DvecType>
pair<shared_ptr<Matrix>, VectorB> FCI_base<CivecType,DvecType>::natorb_convert() {
  assert(rdm1_av_ != nullptr);
  pair<shared_ptr<Matrix>, VectorB> natorb = rdm1_av_->generate_natural_orbitals();
  update_rdms(natorb.first);
  jop_->update_1ext_ints(natorb.first);
  for (auto& i : natorb.second)
    if (i < numerical_zero__) i = 0.0;
  return natorb;
}


template<class CivecType, class DvecType>
void FCI_base<CivecType,DvecType>::update_rdms(shared_ptr<const Matrix> coeff) {
  for (auto& i : *rdm1_)
    i.second->transform(coeff);
  for (auto& i : *rdm2_)
    i.second->transform(coeff);

  // Only when #state > 1, this is needed.
  // Actually rdm1_av_ points to the same object as rdm1_ in 1 state runs. Therefore if you do twice, you get wrong.
  if (rdm1_->size() > 1) rdm1_av_->transform(coeff);
  if (rdm2_->size() > 1) rdm2_av_->transform(coeff);
  assert(rdm1_->size() > 1 || rdm1_->at(0) == rdm1_av_);
  assert(rdm2_->size() > 1 || rdm2_->at(0) == rdm2_av_);
}



template class FCI_base<Civec,Dvec>;
template class FCI_base<DistCivec,DistDvec>;

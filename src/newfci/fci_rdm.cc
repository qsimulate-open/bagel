//
// Newint - Parallel electron correlation program.
// Filename: fci_rdm.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/newfci/fci.h>
#include <src/wfn/rdm.h>

using namespace std;
using namespace bagel;

void NewFCI::compute_rdm12() {
  // Needs initialization here because we use daxpy.
  // For nstate_ == 1, rdm1_av_ = rdm1_[0].
  if (!static_cast<bool>(rdm1_av_) && nstate_ > 1) {
    rdm1_av_ = shared_ptr<RDM<1> >(new RDM<1>(norb_));
    rdm2_av_ = shared_ptr<RDM<2> >(new RDM<2>(norb_));
  }
  if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }
  // we need expanded lists
  shared_ptr<NewDeterminants> detex(new NewDeterminants(norb_, nelea_, neleb_, false));
  cc_->set_det(detex);

  for (int i = 0; i != nstate_; ++i) compute_rdm12(i);

  cc_->set_det(det_);
}


tuple<shared_ptr<RDM<1> >, shared_ptr<RDM<2> > >
  NewFCI::compute_rdm12_last_step(shared_ptr<const NewDvec> dbra, shared_ptr<const NewDvec> dket, shared_ptr<const NewCivec> cibra) const {

  const int nri = dbra->lena()*dbra->lenb();
  const int ij  = norb_*norb_;

  if (nri != dket->lena()*dket->lenb())
    throw logic_error("NewFCI::compute_rdm12_last_step called with inconsistent RI spaces");

  // 1RDM
  // c^dagger <I|\hat{E}|0>
  shared_ptr<RDM<1> > rdm1(new RDM<1>(norb_));
  dgemv_("T", nri, ij, 1.0, dket->data(0)->data(), nri, cibra->data(), 1, 0.0, rdm1->data(), 1);
  // 2RDM
  // \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  shared_ptr<RDM<2> > rdm2(new RDM<2>(norb_));
  dgemm_("T", "N", ij, ij, nri, 1.0, dbra->data(0)->data(), nri, dket->data(0)->data(), nri, 0.0, rdm2->data(), ij);

  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      dcopy_(norb_*norb_, rdm2->element_ptr(0,0,k,i), 1, buf.get(), 1);
      mytranspose1_(buf.get(), &norb_, &norb_, rdm2->element_ptr(0,0,k,i)); // sorting with stride 1 as norb_ is small
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


tuple<shared_ptr<RDM<1> >, shared_ptr<RDM<2> > >
  NewFCI::compute_rdm12_from_civec(shared_ptr<const NewCivec> cbra, shared_ptr<const NewCivec> cket) const {

  // since we consider here number conserving operators...
  shared_ptr<NewDvec> dbra(new NewDvec(cbra->det(), norb_*norb_));
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<NewDvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = shared_ptr<NewDvec>(new NewDvec(cket->det(), norb_*norb_));
    dket->zero();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }

  return compute_rdm12_last_step(dbra, dket, cbra);
}


tuple<shared_ptr<RDM<1> >, shared_ptr<RDM<2> > >
  NewFCI::compute_rdm12_av_from_dvec(shared_ptr<const NewDvec> dbra, shared_ptr<const NewDvec> dket, shared_ptr<const NewDeterminants> o) const {

  if (static_cast<bool>(o)) {
    dbra->set_det(o);
    dket->set_det(o);
  }

  shared_ptr<RDM<1> > rdm1(new RDM<1>(norb_));
  shared_ptr<RDM<2> > rdm2(new RDM<2>(norb_));
  rdm1->zero();
  rdm2->zero();

  assert(dbra->ij() == dket->ij() && dbra->det() == dket->det());

  for (int i = 0; i != dbra->ij(); ++i) {
    shared_ptr<RDM<1> > r1;
    shared_ptr<RDM<2> > r2;
    tie(r1, r2) = compute_rdm12_from_civec(dbra->data(i), dket->data(i));
    rdm1->daxpy(weight_[i], r1);
    rdm2->daxpy(weight_[i], r2);
  }

  if (static_cast<bool>(o)) {
    dbra->set_det(det_);
    dket->set_det(det_);
  }

  return tie(rdm1, rdm2);
}


void NewFCI::compute_rdm12(const int ist) {
  shared_ptr<NewCivec> cc = cc_->data(ist);

  shared_ptr<RDM<1> > rdm1;
  shared_ptr<RDM<2> > rdm2;
  tie(rdm1, rdm2) = compute_rdm12_from_civec(cc, cc);

  // setting to private members.
  rdm1_[ist] = rdm1;
  rdm2_[ist] = rdm2;
  if (nstate_ != 1) {
    rdm1_av_->daxpy(weight_[ist], rdm1);
    rdm2_av_->daxpy(weight_[ist], rdm2);
  } else {
    rdm1_av_ = rdm1;
    rdm2_av_ = rdm2;
  }

}


// note that this does not transform internal integrals (since it is not needed in CASSCF).
pair<vector<double>, vector<double> > NewFCI::natorb_convert() {
  assert(static_cast<bool>(rdm1_av_));
  pair<vector<double>, vector<double> > natorb = rdm1_av_->generate_natural_orbitals();
  update_rdms(natorb.first);
  jop_->update_1ext_ints(natorb.first);
  return natorb;
}


void NewFCI::update_rdms(const vector<double>& coeff) {
  for (auto iter = rdm1_.begin(); iter != rdm1_.end(); ++iter)
    (*iter)->transform(coeff);
  for (auto iter = rdm2_.begin(); iter != rdm2_.end(); ++iter)
    (*iter)->transform(coeff);

  // Only when #state > 1, this is needed.
  // Actually rdm1_av_ points to the same object as rdm1_ in 1 state runs. Therefore if you do twice, you get wrong.
  if (rdm1_.size() > 1) rdm1_av_->transform(coeff);
  if (rdm2_.size() > 1) rdm2_av_->transform(coeff);
}


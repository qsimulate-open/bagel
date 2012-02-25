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


#include <src/fci/fci.h>
#include <src/fci/rdm.h>

using namespace std;

void FCI::compute_rdm12() {
  // Needs initialization here because we use daxpy.
  // For nstate_ == 1, rdm1_av_ = rdm1_[0]. 
  if (!rdm1_av_ && nstate_ > 1) {
    shared_ptr<RDM<1> > rdm1(new RDM<1>(norb_)); 
    shared_ptr<RDM<2> > rdm2(new RDM<2>(norb_)); 
    rdm1_av_ = rdm1;
    rdm2_av_ = rdm2;
    rdm1_av_->zero();
    rdm2_av_->zero();
  }
  for (int i=0; i!=nstate_; ++i) compute_rdm12(i);
}

void FCI::compute_rdm12(const int ist) {
  shared_ptr<Civec> cc = cc_->data(ist);

  // we need expanded lists
  vector<vector<tuple<unsigned int, int, unsigned int> > > phia_bk = phia_;
  vector<vector<tuple<unsigned int, int, unsigned int> > > phib_bk = phib_;
  const_phis_<0>(stringa_, phia_, false);
  const_phis_<1>(stringb_, phib_, false);

  const int la = cc->lena();
  const int lb = cc->lenb();
  const int len = la * lb;
  const int ij = norb_ * norb_;

  // creating new scratch dir.
  shared_ptr<Dvec> d(new Dvec(lb, la, ij));
  d->zero();
  sigma_2a1(cc, d);
  sigma_2a2(cc, d);

  // 1RDM
  shared_ptr<RDM<1> > rdm1(new RDM<1>(norb_));
  dgemv_("T", len, ij, 1.0, d->data(0)->first(), len, cc->first(), 1, 0.0, rdm1->first(), 1);
  // 2RDM
  shared_ptr<RDM<2> > rdm2(new RDM<2>(norb_));
  dgemm_("T", "N", ij, ij, len, 1.0, d->data(0)->first(), len, d->data(0)->first(), len, 0.0, rdm2->first(), ij);

  // sorting... a bit stupid but cheap anyway
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      dcopy_(norb_*norb_, rdm2->element_ptr(0,0,k,i), 1, buf.get(), 1);
      mytranspose1_(buf.get(), &norb_, &norb_, rdm2->element_ptr(0,0,k,i)); // sorting with stride 1 as norb_ is small
    }
  }

  // put int diagonal into 2RDM
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      for (int j = 0; j != norb_; ++j) {
        rdm2->element(j,k,k,i) -= rdm1->element(j,i); 
      }
    }
  }

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

  phia_ = phia_bk;
  phib_ = phib_bk;
}


// note that this does not transform internal integrals (since it is not needed in CASSCF). 
pair<vector<double>, vector<double> > FCI::natorb_convert() {
  assert(rdm1_av_);
  pair<vector<double>, vector<double> > natorb = rdm1_av_->generate_natural_orbitals(); 
  update_rdms(natorb.first);
  jop_->update_1ext_ints(natorb.first);
  return natorb;
}


void FCI::update_rdms(const vector<double>& coeff) {
  for (auto iter = rdm1_.begin(); iter != rdm1_.end(); ++iter)
    (*iter)->transform(coeff);
  for (auto iter = rdm2_.begin(); iter != rdm2_.end(); ++iter)
    (*iter)->transform(coeff);

  // Only when #state > 1, this is needed.
  // Actually rdm1_av_ points to the same object as rdm1_ in 1 state runs. Therefore if you do twice, you get wrong.
  if (rdm1_.size() > 1) rdm1_av_->transform(coeff);
  if (rdm2_.size() > 1) rdm2_av_->transform(coeff);
}


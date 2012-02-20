//
// Newint - Parallel electron correlation program.
// Filename: superci_sigma.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#include <src/casscf/superci.h>

using namespace std;


//////////////////////////////////////////////////////
// gradient vectors ... they are definitely correct.
//////////////////////////////////////////////////////


// <a/i|H|0> = 2f_ai /sqrt(2)
void SuperCI::grad_vc(const shared_ptr<Matrix1e> f, shared_ptr<RotFile> sigma) {
  if (!nvirt_ || !nclosed_) return;
  double* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_; ++i, target += nvirt_)
    daxpy_(nvirt_, std::sqrt(2.0), f->element_ptr(nocc_,i), 1, target, 1);
}


// <a/r|H|0> finact_as d_sr + 2(as|tu)P_rs,tu = fact_ar  (/sqrt(nr))
void SuperCI::grad_va(const shared_ptr<QFile> fact, shared_ptr<RotFile> sigma) {
  if (!nvirt_ || !nact_) return;
  double* target = sigma->ptr_va();
  for (int i = 0; i != nact_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 1.0/std::sqrt(occup_[i]), fact->element_ptr(nocc_, i), 1, target, 1);
  }
}


// <r/i|H|0> = (2f_ri - f^act_ri)/sqrt(2-nr)
void SuperCI::grad_ca(const shared_ptr<Matrix1e> f, shared_ptr<QFile> fact, shared_ptr<RotFile> sigma) {
  if (!nclosed_ || !nact_) return;
  double* target = sigma->ptr_ca();
  for (int i = 0; i != nact_; ++i, target += nclosed_) {
    daxpy_(nclosed_, 2.0/std::sqrt(2.0-occup_[i]), f->element_ptr(0,nclosed_+i), 1, target, 1);
    daxpy_(nclosed_, -1.0/std::sqrt(2.0-occup_[i]), fact->element_ptr(0,i), 1, target, 1);
  }
}



//////////////////////////////////////////////////////
// other elements ... hopefully correct.
//////////////////////////////////////////////////////


// sigma_at_at = delta_ab Gtu/sqrt(nt nu) + delta_tu Fab
void SuperCI::sigma_at_at_(const shared_ptr<RotFile> cc, shared_ptr<RotFile> sigma, const shared_ptr<QFile> gaa, const shared_ptr<Matrix1e> f) {
  if (!nact_ || !nvirt_) return;
  shared_ptr<QFile> gtup(new QFile(*gaa));
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      gtup->element(j,i) /= std::sqrt(occup_[i]*occup_[j]);
  dgemm_("N", "N", nvirt_, nact_, nact_, 1.0, cc->ptr_va(), nvirt_, gtup->data(), nact_, 1.0, sigma->ptr_va(), nvirt_);
  dgemm_("N", "N", nvirt_, nact_, nvirt_, 1.0, f->element_ptr(nocc_, nocc_), nbasis_, cc->ptr_va(), nvirt_, 1.0, sigma->ptr_va(), nvirt_);
}


// sigma_ai_ai = delta_ij F_ab - delta_ab F_ij
void SuperCI::sigma_ai_ai_(const shared_ptr<RotFile> cc, shared_ptr<RotFile> sigma, const shared_ptr<Matrix1e> f) {
  if (!nact_ || !nvirt_) return;
  dgemm_("N", "N", nvirt_, nclosed_, nclosed_, -1.0, cc->ptr_vc(), nvirt_, f->data(), nbasis_, 1.0, sigma->ptr_vc(), nvirt_);
  dgemm_("N", "N", nvirt_, nclosed_, nvirt_, 1.0, f->element_ptr(nocc_, nocc_), nbasis_, cc->ptr_vc(), nvirt_, 1.0, sigma->ptr_vc(), nvirt_);
}


// sigma_at_ai = -delta_ab Fact_ti sqrt(nt/2)
void SuperCI::sigma_at_ai_(const shared_ptr<RotFile> cc, shared_ptr<RotFile> sigma, const shared_ptr<QFile> fact) {
  if (!nact_ || !nvirt_ || !nclosed_) return;
  QFile tmp(nclosed_, nact_);
  tmp.zero();
  for (int i = 0; i != nact_; ++i) {
    const double fac = -std::sqrt(0.5*occup_[i]);
    daxpy_(nclosed_, fac, fact->element_ptr(0,i), 1, tmp.element_ptr(0,i), 1); 
  }
  dgemm_("N", "N", nvirt_, nact_, nclosed_, 1.0, cc->ptr_vc(), nvirt_, tmp.data(), nclosed_, 1.0, sigma->ptr_va(), nvirt_); 
  dgemm_("N", "T", nvirt_, nclosed_, nact_, 1.0, cc->ptr_va(), nvirt_, tmp.data(), nclosed_, 1.0, sigma->ptr_vc(), nvirt_); 
}


// sigma_ai_ti = sqrt((2-nt)/2) Fact_at
void SuperCI::sigma_ai_ti_(const shared_ptr<RotFile> cc, shared_ptr<RotFile> sigma, const shared_ptr<QFile> fact) {
  if (!nact_ || !nvirt_ || !nclosed_) return;
  QFile tmp(nvirt_, nact_);
  tmp.zero();
  for (int i = 0; i != nact_; ++i) {
    const double fac = std::sqrt(1.0-0.5*occup_[i]);
    daxpy_(nvirt_, fac, fact->element_ptr(nocc_,i), 1, tmp.element_ptr(0,i), 1); 
  }
  dgemm_("T", "N", nclosed_, nact_, nvirt_, 1.0, cc->ptr_vc(), nvirt_, tmp.data(), nvirt_, 1.0, sigma->ptr_ca(), nclosed_); 
  dgemm_("N", "T", nvirt_, nclosed_, nact_, 1.0, tmp.data(), nvirt_, cc->ptr_ca(), nclosed_, 1.0, sigma->ptr_vc(), nvirt_);

}


// sigma_ti_ti = - delta_ij ((2-nt-nu)Fact_tu - G_tu)/sqrt((2-nt)(2-nu)) - delta_tu f_ij
void SuperCI::sigma_ti_ti_(const shared_ptr<RotFile> cc, shared_ptr<RotFile> sigma, const shared_ptr<QFile> gaa, const shared_ptr<Matrix1e> f,
                           const shared_ptr<QFile> factp) {
  if (!nact_ || !nclosed_) return;
  QFile tmp(nact_, nact_);
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nact_; ++j) {
      tmp.element(j,i) = -((2.0 - occup_[j] - occup_[i]) * factp->element(j,i) - gaa->element(j,i))
                         / std::sqrt((2.0 - occup_[i]) * (2.0 - occup_[j]));
    }
  }
  dgemm_("N", "N", nclosed_, nact_, nact_, 1.0, cc->ptr_ca(), nclosed_, tmp.data(), nact_, 1.0, sigma->ptr_ca(), nclosed_); 
  dgemm_("N", "N", nclosed_, nact_, nclosed_, -1.0, f->data(), nbasis_, cc->ptr_ca(), nclosed_, 1.0, sigma->ptr_ca(), nclosed_); 
}


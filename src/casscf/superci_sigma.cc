//
// BAGEL - Parallel electron correlation program.
// Filename: superci_sigma.cc
// Copyright (C) 2011 Toru Shiozaki
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


#include <src/casscf/superci.h>

using namespace std;
using namespace bagel;


//////////////////////////////////////////////////////
// gradient vectors ... they are definitely correct.
//////////////////////////////////////////////////////


// <a/i|H|0> = 2f_ai /sqrt(2)
void SuperCI::grad_vc(const shared_ptr<Matrix> f, shared_ptr<RotFile> sigma) {
  if (!nvirt_ || !nclosed_) return;
  double* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_; ++i, target += nvirt_)
    daxpy_(nvirt_, std::sqrt(2.0), f->element_ptr(nocc_,i), 1, target, 1);
}


// <a/r|H|0> finact_as h_sr + (as|tu)D_rs,tu = fact_ar  (/sqrt(nr) - due to normalization)
void SuperCI::grad_va(const shared_ptr<Matrix> fact, shared_ptr<RotFile> sigma) {
  if (!nvirt_ || !nact_) return;
  double* target = sigma->ptr_va();
  for (int i = 0; i != nact_; ++i, target += nvirt_) {
    const double fac = (occup_[i]>occup_thresh) ? 1.0/std::sqrt(occup_[i]) : 0.0;
    daxpy_(nvirt_, fac, fact->element_ptr(nocc_, i), 1, target, 1);
  }
}


// <r/i|H|0> = (2f_ri - f^act_ri)/sqrt(2-nr)
void SuperCI::grad_ca(const shared_ptr<Matrix> f, shared_ptr<Matrix> fact, shared_ptr<RotFile> sigma) {
  if (!nclosed_ || !nact_) return;
  double* target = sigma->ptr_ca();
  for (int i = 0; i != nact_; ++i, target += nclosed_) {
    const double fac = (2.0-occup_[i] > occup_thresh) ? 1.0/std::sqrt(2.0-occup_[i]) : 0.0;
    daxpy_(nclosed_, 2.0*fac, f->element_ptr(0,nclosed_+i), 1, target, 1);
    daxpy_(nclosed_, -fac, fact->element_ptr(0,i), 1, target, 1);
  }
}




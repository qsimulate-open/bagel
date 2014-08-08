//
// BAGEL - Parallel electron correlation program.
// Filename: zsuperci_sigma.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson Bates <jefferson.bates@northwestern.edu>
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


#include <src/zcasscf/zsuperci.h>

using namespace std;
using namespace bagel;

// <a/i|H|0> = f_ai
void ZSuperCI::grad_vc(const shared_ptr<ZMatrix> f, shared_ptr<ZRotFile> sigma) {
#ifdef BOTHSPACES
  const int nvirt_tmp = nvirt_;
#else
  const int nvirt_tmp = nvirtnr_;
#endif
  if (!nvirt_tmp || !nclosed_) return;
  complex<double>* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_*2; ++i, target += nvirt_tmp*2)
    zaxpy_(nvirt_tmp*2, 1.0, f->element_ptr(nocc_*2,i), 1, target, 1);
}


// <a/r|H|0> = cfock_as^* n_s + (as|tu)D_rs,tu = fact_ar  (/sqrt(n_s) - due to normalization)
void ZSuperCI::grad_va(const shared_ptr<ZMatrix> fact, shared_ptr<ZRotFile> sigma) {
#ifdef BOTHSPACES
  const int nvirt_tmp = nvirt_;
#else
  const int nvirt_tmp = nvirtnr_;
#endif
  if (!nvirt_ || !nact_) return;
  complex<double>* target = sigma->ptr_va();
  for (int i = 0; i != nact_*2; ++i, target += nvirt_tmp*2) {
    const double fac = (occup_[i]>zoccup_thresh) ? 1.0/std::sqrt(occup_[i]) : 0.0;
    zaxpy_(nvirt_tmp*2, fac, fact->element_ptr(nocc_*2, i), 1, target, 1);
  }
}


// <r/i|H|0> = (f_ri - f^act_ri)/sqrt(1-nr)
void ZSuperCI::grad_ca(const shared_ptr<ZMatrix> f, shared_ptr<ZMatrix> fact, shared_ptr<ZRotFile> sigma) {
  if (!nclosed_ || !nact_) return;
  complex<double>* target = sigma->ptr_ca();
  shared_ptr<ZMatrix> f_conj = f->slice_copy(nclosed_*2, nocc_*2)->get_conjg();
  for (int i = 0; i != nact_*2; ++i, target += nclosed_*2) {
    const double fac = (1.0-occup_[i] > zoccup_thresh) ? 1.0/std::sqrt(1.0-occup_[i]) : 0.0;
    zaxpy_(nclosed_*2, fac, f_conj->element_ptr(0,i), 1, target, 1);
    zaxpy_(nclosed_*2, -fac, fact->element_ptr(0,i), 1, target, 1);
  }
}


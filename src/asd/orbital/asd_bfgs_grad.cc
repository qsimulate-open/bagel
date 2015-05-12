//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/asd_bfgs_grad.hpp
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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

#include <src/asd/orbital/asd_bfgs.h>

using namespace std;
using namespace bagel;

// grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
// G_(i<a)
void ASD_BFGS::grad_vc(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<ASD_RotFile> grad) const {
  if (!nvirt_ || !nclosed_) return;
  if (!grad->is_inter() && grad->is_intra()) return;
  double* target = grad->ptr_vc();
  for (int i = 0; i != nclosed_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 4.0, cfock->element_ptr(nocc_,i), 1, target, 1);
    daxpy_(nvirt_, 4.0, afock->element_ptr(nocc_,i), 1, target, 1);
  }
}

// grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
// G_(t<a)
void ASD_BFGS::grad_va(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> qxr, shared_ptr<Matrix> rdm1, shared_ptr<ASD_RotFile> grad) const {
  if (!nvirt_ || !nact_) return;
  if (!grad->is_inter() && grad->is_intra()) return;
  dgemm_("N", "T", nvirt_, nact_, nact_, 2.0, cfock->element_ptr(nocc_,nclosed_), cfock->ndim(), rdm1->data(), rdm1->ndim(), 0.0, grad->ptr_va(), nvirt_);
  double* target = grad->ptr_va();
  for (int i = 0; i != nact_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 2.0, qxr->element_ptr(nocc_, i), 1, target, 1);
  }
}

// grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
// G_(i<t)
void ASD_BFGS::grad_ca(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<Matrix> rdm1, shared_ptr<ASD_RotFile> grad) const {
  if (!nclosed_ || !nact_) return;
  if (!grad->is_inter() && grad->is_intra()) return;
  double* target = grad->ptr_ca();
  for (int i = 0; i != nact_; ++i, target += nclosed_) {
    daxpy_(nclosed_, 4.0, cfock->element_ptr(0,nclosed_+i), 1, target, 1);
    daxpy_(nclosed_, 4.0, afock->element_ptr(0,nclosed_+i), 1, target, 1);
    daxpy_(nclosed_, -2.0, qxr->element_ptr(0, i), 1, target, 1);
  }
  //-2 cfock_iu * D_ur
  dgemm_("T", "N", nclosed_, nact_, nact_, -2.0, cfock->element_ptr(nclosed_,0), cfock->ndim(), rdm1->data(), rdm1->ndim(), 1.0, grad->ptr_ca(), nclosed_);
}

// grad(t/t)
// G_(t(A)<t(B))
void ASD_BFGS::grad_aa(shared_ptr<const Matrix> mcfock, shared_ptr<ASD_RotFile> grad) const {
  if (grad->is_inter() && !grad->is_intra()) return;
  if (!nact_) return;
  {
    double* target = grad->ptr_aa();
    for (int ia = 0; ia != nactA_; ++ia)  //A
      for (int jb = nactA_; jb != nact_; ++jb, ++target)  //B
        *target = 2.0*(mcfock->element(jb,ia) - mcfock->element(ia,jb));
  }
  //RAS
  if (nactA_) { //Monomer A
    if (rasA_[0]) {//21
      double* target = grad->ptr_aa21A();
      for (int ia = 0; ia != rasA_[0]; ++ia)  //1
        for (int jb = rasA_[0]; jb != rasA_[0]+rasA_[1]; ++jb, ++target)  //2
          *target = 2.0*(mcfock->element(jb,ia) - mcfock->element(ia,jb));
    }
    if (rasA_[0] && rasA_[2]) {//31
      double* target = grad->ptr_aa31A();
      for (int ia = 0; ia != rasA_[0]; ++ia)  //1
        for (int jb = rasA_[0]+rasA_[1]; jb != nactA_; ++jb, ++target)  //3
          *target = 2.0*(mcfock->element(jb,ia) - mcfock->element(ia,jb));
    }
    if (rasA_[2]) {//32
      double* target = grad->ptr_aa32A();
      for (int ia = rasA_[0]; ia != rasA_[0]+rasA_[1]; ++ia) //2
        for (int jb = rasA_[0]+rasA_[1]; jb != nactA_; ++jb, ++target) //3
          *target = 2.0*(mcfock->element(jb,ia) - mcfock->element(ia,jb));
    }
  }

  if (nactB_) {//Monomer B :21
    if (rasB_[0]) {
      double* target = grad->ptr_aa21B();
      for (int ia = nactA_; ia != nactA_+rasB_[0]; ++ia)  //1
        for (int jb = nactA_+rasB_[0]; jb != nactA_+rasB_[0]+rasB_[1]; ++jb, ++target)  //2
          *target = 2.0*(mcfock->element(jb,ia) - mcfock->element(ia,jb));
    }
    if (rasB_[0] && rasB_[2]) {//31
      double* target = grad->ptr_aa31B();
      for (int ia = nactA_; ia != nactA_+rasB_[0]; ++ia)  //1
        for (int jb = nactA_+rasB_[0]+rasB_[1]; jb != nact_; ++jb, ++target)  //3
          *target = 2.0*(mcfock->element(jb,ia) - mcfock->element(ia,jb));
    }
    if (rasB_[2]) {//32
      double* target = grad->ptr_aa32B();
      for (int ia = nactA_+rasB_[0]; ia != nactA_+rasB_[0]+rasB_[1]; ++ia)  //2
        for (int jb = nactA_+rasB_[0]+rasB_[1]; jb != nact_; ++jb, ++target)  //3
          *target = 2.0*(mcfock->element(jb,ia) - mcfock->element(ia,jb));
    }
  }
}


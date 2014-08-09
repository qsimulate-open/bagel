//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf_denom.cc
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

#include <src/zcasscf/zcasbfgs.h>

using namespace std;
using namespace bagel;

shared_ptr<ZRotFile> ZCASBFGS::compute_denom(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1) const {
  auto out = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);

  shared_ptr<ZMatrix> cfockd;
  if (nact_) {
    cfockd = make_shared<ZMatrix>(*cfock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2) * *rdm1);
    cfockd->hermite();
  }

  // ia part (4.7a)
  if (nvirt_ && nclosed_) {
    complex<double>* target = out->ptr_vc();
    for (int i = 0; i != nclosed_*2; ++i) {
      for (int j = 0; j != nvirt_*2; ++j) {
        *target++ = cfock->element(j+nocc_*2, j+nocc_*2) + afock->element(j+nocc_*2, j+nocc_*2) - cfock->element(i,i) - afock->element(i,i);
      }
    }
  }
  // ra part (4.7b)
  if (nvirt_ && nact_) {
    complex<double>* target = out->ptr_va();
    for (int i = 0; i != nact_*2; ++i) {
      for (int j = 0; j != nvirt_*2; ++j) {
        *target++ = rdm1->element(i, i)*(cfock->element(j+nocc_*2, j+nocc_*2)+afock->element(j+nocc_*2, j+nocc_*2))
                  - cfockd->element(i, i) - qxr->element(i+nclosed_*2, i);
      }
    }
  }
  // it part (4.7c)
  if (nclosed_ && nact_) {
    complex<double>* target = out->ptr_ca();
    for (int i = 0; i != nact_*2; ++i) {
      for (int j = 0; j != nclosed_*2; ++j) {
        *target++ = (cfock->element(i+nclosed_*2, i+nclosed_*2)+afock->element(i+nclosed_*2, i+nclosed_*2) - cfock->element(j,j) - afock->element(j,j))
                  + rdm1->element(i,i)*(cfock->element(j,j)+afock->element(j,j)) - cfockd->element(i,i) - qxr->element(i+nclosed_*2, i);
      }
    }
  }

  const double thresh = 1.0e-8;
  for (int i = 0; i != out->size(); ++i)
    if (fabs(out->data(i)) < thresh) {
      out->data(i) = 1.0e10;
    }
  return out;
}



// grad(a/i) (eq.4.3a): (cfock_ai+afock_ai)
void ZCASBFGS::grad_vc(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<ZRotFile> sigma) const {
  if (!nvirt_ || !nclosed_) return;
  complex<double>* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_*2; ++i, target += nvirt_*2) {
    zaxpy_(nvirt_*2, 1.0, cfock->element_ptr(nocc_*2, i), 1, target, 1);
    zaxpy_(nvirt_*2, 1.0, afock->element_ptr(nocc_*2, i), 1, target, 1);
  }
}


// grad(a/t) (eq.4.3b): cfock_au gamma_ut + q_at^*
void ZCASBFGS::grad_va(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1, shared_ptr<ZRotFile> sigma) const {
  if (!nvirt_ || !nact_) return;
  zgemm3m_("N", "N", nvirt_*2, nact_*2, nact_*2, 1.0, cfock->element_ptr(nocc_*2, nclosed_*2), cfock->ndim(), rdm1->data(), rdm1->ndim(), 0.0, sigma->ptr_va(), nvirt_*2);
  complex<double>* target = sigma->ptr_va();
  for (int i = 0; i != nact_*2; ++i, target += nvirt_*2) {
    zaxpy_(nvirt_*2, 1.0, qxr->element_ptr(nocc_*2, i), 1, target, 1);
  }
}


// grad(r/i) (eq.4.3c): (cfock_ri+afock_ri) - cfock_iu gamma_ur - qxr_ir
void ZCASBFGS::grad_ca(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1, shared_ptr<ZRotFile> sigma) const {
  if (!nclosed_ || !nact_) return;
  auto qxrc = qxr->get_conjg();
  auto afockc = afock->get_conjg();
  auto cfockc = cfock->get_conjg();
  complex<double>* target = sigma->ptr_ca();
  for (int i = 0; i != nact_*2; ++i, target += nclosed_*2) {
    zaxpy_(nclosed_*2, 1.0, afockc->element_ptr(0,nclosed_*2+i), 1, target, 1);
    zaxpy_(nclosed_*2, 1.0, cfockc->element_ptr(0,nclosed_*2+i), 1, target, 1);
    zaxpy_(nclosed_*2, -1.0, qxrc->element_ptr(0, i), 1, target, 1);
  }
  // "T" effectively makes complex conjugate of cfock
  zgemm3m_("T", "N", nclosed_*2, nact_*2, nact_*2, -1.0, cfock->element_ptr(nclosed_*2, 0), cfock->ndim(), rdm1->data(), rdm1->ndim(), 1.0, sigma->ptr_ca(), nclosed_*2);
}

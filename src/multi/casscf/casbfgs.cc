//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: casbfgs_base.cc
// Copyright (C) 2013 Toru Shiozaki
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


#include <src/multi/casscf/casbfgs.h>

using namespace std;
using namespace bagel;


void CASBFGS::compute() {
  const string algo = idata_->get<string>("bfgstype", "bagel");
  if (algo != "bagel" && algo != "alglib")
    throw runtime_error("unknown BFGS type specified");
  // first do BFGS2 (alglib version)
  if (algo != "bagel") {
    auto bfgs = make_shared<CASBFGS2>(idata_, geom_, ref_);
    bfgs->compute();
    refout_ = bfgs->conv_to_ref();

    if (!bfgs->only_energy_converged()) {
      fci_ = bfgs->fci();
      energy_ = bfgs->energy();
      rms_grad_ = bfgs->rms_grad();
      return;
    }
  }
  // second do BAGEL's native BFGS
  {
    auto bfgs = make_shared<CASBFGS1>(idata_, geom_, refout_);
    bfgs->compute();
    refout_ = bfgs->conv_to_ref();

    fci_ = bfgs->fci();
    energy_ = bfgs->energy();
    rms_grad_ = bfgs->rms_grad();
  }
}


shared_ptr<const RotFile> CASBFGS_base::compute_denom(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {
  auto out = make_shared<RotFile>(nclosed_, nact_, nvirt_);
  const double tiny = 1.0e-15;

  // ia part (4.7a)
  if (nvirt_ && nclosed_) {
    double* target = out->ptr_vc();
    for (int i = 0; i != nclosed_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        *target++ = 4.0*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_)) - 4.0*(cfock->element(i,i)+afock->element(i,i));
      }
    }
  }
  // ra part (4.7b)
  if (nvirt_ && nact_) {
    double* target = out->ptr_va();
    for (int i = 0; i != nact_; ++i) {
      if (occup_[i] < tiny) continue;
      for (int j = 0; j != nvirt_; ++j) {
        *target++ = 2.0*occup_[i]*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_))
                  - 2.0*occup_[i]*cfock->element(i+nclosed_, i+nclosed_) - 2.0*qxr->element(i+nclosed_, i);
      }
    }
  }
  // it part (4.7c)
  if (nclosed_ && nact_) {
    double* target = out->ptr_ca();
    for (int i = 0; i != nact_; ++i) {
      if (occup_[i] < tiny) continue;
      for (int j = 0; j != nclosed_; ++j) {
        *target++ = 4.0*(cfock->element(i+nclosed_, i+nclosed_)+afock->element(i+nclosed_, i+nclosed_) - cfock->element(j,j) - afock->element(j,j))
                  + 2.0*occup_[i]*(cfock->element(j,j)+afock->element(j,j)) - 2.0*occup_[i]*cfock->element(i+nclosed_, i+nclosed_) - 2.0*qxr->element(i+nclosed_, i);
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


// grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
void CASBFGS_base::grad_vc(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<RotFile> sigma) const {
  if (!nvirt_ || !nclosed_) return;
  double* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 4.0, cfock->element_ptr(nocc_,i), 1, target, 1);
    daxpy_(nvirt_, 4.0, afock->element_ptr(nocc_,i), 1, target, 1);
  }
}


// grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
// gamma is assumed to be diagonal
void CASBFGS_base::grad_va(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> qxr, shared_ptr<RotFile> sigma) const {
  if (!nvirt_ || !nact_) return;
  shared_ptr<const RDM<1>> rdm1 = fci_->rdm1_av();
  double* target = sigma->ptr_va();
  for (int i = 0; i != nact_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 2.0, qxr->element_ptr(nocc_, i), 1, target, 1);
    for (int j = 0; j != nact_; ++j)
      daxpy_(nvirt_, 2.0*rdm1->element(j,i), cfock->element_ptr(nocc_, nclosed_+j), 1, target, 1);
  }
}


// grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
void CASBFGS_base::grad_ca(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<RotFile> sigma) const {
  if (!nclosed_ || !nact_) return;
  shared_ptr<const RDM<1>> rdm1 = fci_->rdm1_av();
  {
    double* target = sigma->ptr_ca();
    for (int i = 0; i != nact_; ++i, target += nclosed_) {
      daxpy_(nclosed_, 4.0, cfock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, 4.0, afock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, -2.0, qxr->element_ptr(0, i), 1, target, 1);
      for (int j = 0; j != nact_; ++j)
        daxpy_(nclosed_, -2.0*rdm1->element(j,i), cfock->element_ptr(0,nclosed_+j), 1, target, 1);
    }
  }
}


//
// BAGEL - Parallel electron correlation program.
// Filename: casbfgs.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <src/casscf/casbfgs.h>
#include <src/casscf/qvec.h>
#include <iostream>
#include <src/fci/fci.h>
#include <src/util/davidson.h>
#include <src/scf/hcore.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <src/util/bfgs.h>

using namespace std;
using namespace bagel;

void CASBFGS::compute() {

  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.

  shared_ptr<BFGS<RotFile>> bfgs;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  shared_ptr<Matrix> x(new Matrix(nbasis_, nbasis_));
  x->unit();

  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    mute_stdcout();
    if (iter) fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    energy_ = fci_->energy();
    resume_stdcout();

    // here make a natural orbitals and update the coefficients. Closed and virtual orbitals remain canonical
    shared_ptr<Matrix> natorb = form_natural_orbs();

    shared_ptr<RotFile> sigma_(new RotFile(nclosed_, nact_, nvirt_, false));
    sigma_->zero();

    // compute one-boedy operators
    // * preparation
    shared_ptr<const Matrix> ccoeff = coeff_->slice(0, nclosed_);
    shared_ptr<const Matrix> ocoeff = coeff_->slice(0, nocc_);
    // * core Fock operator
    shared_ptr<const Matrix> cden = nclosed_ ? coeff_->form_density_rhf(nclosed_, 0) : shared_ptr<const Matrix>(new Matrix(nbasis_, nbasis_));
    shared_ptr<const Matrix> cfockao = nclosed_ ? shared_ptr<const Matrix>(new Fock<1>(geom_, hcore_, cden, ccoeff)) : hcore_;
    shared_ptr<const Matrix> cfock(new Matrix(*coeff_ % *cfockao * *coeff_));
    // * active Fock operator
    // first make a weighted coefficient
    shared_ptr<Matrix> acoeff = coeff_->slice(nclosed_, nocc_);
    for (int i = 0; i != nact_; ++i)
      dscal_(acoeff->ndim(), sqrt(occup_[i]/2.0), acoeff->element_ptr(0, i), 1); 
    // then make a AO density matrix
    shared_ptr<const Matrix> aden(new Matrix((*acoeff ^ *acoeff)*2.0));
    shared_ptr<const Matrix> afockao(new Fock<1>(geom_, hcore_, aden, acoeff));
    shared_ptr<const Matrix> afock(new Matrix(*coeff_ % (*afockao - *hcore_) * *coeff_));

    // * Q_xr = 2(xs|tu)P_rs,tu (x=general, mo)
    shared_ptr<const Matrix> qxr(new Qvec(geom_->nbasis(), nact_, geom_->df(), coeff_, nclosed_, fci_, fci_->rdm2_av()));

    // grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
    grad_vc(cfock, afock, sigma_);
    // grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at 
    grad_va(cfock, qxr, sigma_);
    // grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir 
    grad_ca(cfock, afock, qxr, sigma_);

    // if this is the first time, set up the BFGS solver
    if (iter == 0) {
      shared_ptr<const RotFile> denom = compute_denom(cfock, afock, qxr);
      bfgs = shared_ptr<BFGS<RotFile>>(new BFGS<RotFile>(denom));
    }
    // extrapolation using BFGS
    shared_ptr<RotFile> xrot(new RotFile(x->log(), nclosed_, nact_, nvirt_, false));
    shared_ptr<RotFile> a = bfgs->extrapolate(sigma_, xrot);
    *a *= -1.0;

    // restore the matrix from RotFile
    shared_ptr<const Matrix> amat = a->unpack();
    shared_ptr<Matrix> expa = amat->exp(2);
    expa->purify_unitary();
    coeff_ = shared_ptr<const Coeff>(new Coeff(*coeff_**expa));

    // for next BFGS extrapolation
    *x *= *expa;

    // setting error of macro iteration
    const double gradient = sigma_->ddot(*sigma_) / sigma_->size();

    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    if (gradient < thresh_) break;

    if (iter == max_iter_-1)
      throw runtime_error("Max iteration reached in the CASSCF macro interation.");
  }
  // ============================
  // macro iteration to here
  // ============================

#if 0
  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  fci_->update(coeff_);
  fci_->compute();
  fci_->compute_rdm12();
#endif

  throw logic_error("not yet fully implemented");

}


shared_ptr<const RotFile> CASBFGS::compute_denom(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {
  shared_ptr<RotFile> out(new RotFile(nclosed_, nact_, nvirt_, false));
  out->fill(1.0);
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
    if (fabs(out->data(i)) < thresh)
      out->data(i) = 1.0e10;
  return out;
}


// grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
void CASBFGS::grad_vc(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<RotFile> sigma) const {
  if (!nvirt_ || !nclosed_) return;
  double* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 4.0, cfock->element_ptr(nocc_,i), 1, target, 1);
    daxpy_(nvirt_, 4.0, afock->element_ptr(nocc_,i), 1, target, 1);
  }
}


// grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at 
// gamma is assumed to be diagonal
void CASBFGS::grad_va(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> qxr, shared_ptr<RotFile> sigma) const {
  if (!nvirt_ || !nact_) return;
  double* target = sigma->ptr_va();
  for (int i = 0; i != nact_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 2.0*occup_[i], cfock->element_ptr(nocc_, i+nclosed_), 1, target, 1);
    daxpy_(nvirt_, 2.0, qxr->element_ptr(nocc_, i), 1, target, 1);
  }
}


// grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir 
void CASBFGS::grad_ca(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<RotFile> sigma) const {
  if (!nclosed_ || !nact_) return;
  {
    double* target = sigma->ptr_ca();
    for (int i = 0; i != nact_; ++i, target += nclosed_) {
      daxpy_(nclosed_, 4.0-2.0*occup_[i], cfock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, 4.0, afock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, -2.0, qxr->element_ptr(0, i), 1, target, 1); 
    }
  }
}


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


#include <src/casscf/casbfgs.h>
#include <src/casscf/qvec.h>
#include <src/math/davidson.h>
#include <src/math/bfgs.h>
#include <src/math/hpw_diis.h>

using namespace std;
using namespace bagel;

void CASBFGS::compute() {

  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.

  shared_ptr<BFGS<Matrix>> bfgs;
  shared_ptr<HPW_DIIS<Matrix>> diis;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  auto x = make_shared<Matrix>(nbasis_, nbasis_);
  x->unit();
  shared_ptr<const Matrix> xstart;

  for (int iter = 0; iter != max_iter_; ++iter) {

    const shared_ptr<const Coeff> cold = coeff_;
    const shared_ptr<const Matrix> xold = x->copy();

    // first perform CASCI to obtain RDMs
    mute_stdcout();
    if (iter) fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    energy_ = fci_->energy();
    resume_stdcout();

    shared_ptr<Matrix> natorb_mat = x->clone();
    {
      // here make a natural orbitals and update coeff_. Closed and virtual orbitals remain canonical. Also, FCI::rdms are updated
      shared_ptr<const Matrix> natorb = form_natural_orbs();
      natorb_mat->unit();
      natorb_mat->copy_block(nclosed_, nclosed_, nact_, nact_, natorb);
    }

    auto sigma = make_shared<RotFile>(nclosed_, nact_, nvirt_, false);
    sigma->zero();

    // compute one-boedy operators
    // * preparation
    shared_ptr<const Matrix> ccoeff = coeff_->slice(0, nclosed_);
    shared_ptr<const Matrix> ocoeff = coeff_->slice(0, nocc_);
    // * core Fock operator
    shared_ptr<const Matrix> cden = nclosed_ ? coeff_->form_density_rhf(nclosed_, 0) : make_shared<const Matrix>(geom_->nbasis(), geom_->nbasis());
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<const Fock<1>>(geom_, hcore_, cden, ccoeff) : hcore_;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    // * active Fock operator
    // first make a weighted coefficient
    shared_ptr<Matrix> acoeff = coeff_->slice(nclosed_, nocc_);
    for (int i = 0; i != nact_; ++i)
      dscal_(acoeff->ndim(), sqrt(occup_[i]/2.0), acoeff->element_ptr(0, i), 1);
    // then make a AO density matrix
    shared_ptr<const Matrix> aden = make_shared<Matrix>((*acoeff ^ *acoeff)*2.0);
    shared_ptr<const Matrix> afockao = make_shared<Fock<1>>(geom_, hcore_, aden, acoeff);
    shared_ptr<const Matrix> afock = make_shared<Matrix>(*coeff_ % (*afockao - *hcore_) * *coeff_);

    // * Q_xr = 2(xs|tu)P_rs,tu (x=general, mo)
    auto qxr = make_shared<const Qvec>(coeff_->mdim(), nact_, geom_->df(), coeff_, nclosed_, fci_, fci_->rdm2_av());

    // grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
    grad_vc(cfock, afock, sigma);
    // grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
    grad_va(cfock, qxr, sigma);
    // grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
    grad_ca(cfock, afock, qxr, sigma);

    // if this is the first time, set up the BFGS solver
//  if (iter == 0) {
  if (true) {
      // BFGS and DIIS should start at the same time
      shared_ptr<const Matrix> denom = compute_denom(cfock, afock, qxr)->unpack(1.0e10);
      bfgs = make_shared<BFGS<Matrix>>(denom);
    }
    if (iter == 0) {
//if (false) {
      xstart = xold->copy();
      shared_ptr<Matrix> unit = make_shared<Matrix>(xold->mdim(), xold->mdim());
      unit->unit();
      diis = make_shared<HPW_DIIS<Matrix>>(10, cold, unit);
    }
    // extrapolation using BFGS
    *x *= *natorb_mat;
    auto xlog = make_shared<Matrix>(*x->log(100));
    shared_ptr<const Matrix> sigma_mat = sigma->unpack();
    shared_ptr<Matrix> a = bfgs->extrapolate(sigma_mat, xlog);
    *a *= -1.0;

    // restore the matrix from RotFile
    shared_ptr<const Matrix> amat = a;
    shared_ptr<Matrix> expa = amat->exp(100);
    expa->purify_unitary();

    if (!diis) {
      coeff_ = make_shared<const Coeff>(*coeff_**expa);
      // for next BFGS extrapolation
      *x *= *expa;
    } else {
      auto tmp3 = make_shared<const Matrix>(*natorb_mat * *expa ^ *natorb_mat);
      shared_ptr<const Matrix> mcc = diis->extrapolate(tmp3);
      coeff_ = make_shared<const Coeff>(*mcc);
      // update x
      x = make_shared<Matrix>(*xstart * *diis->extrap());
//    cout << setprecision(10) << (*coeff_ - *diis->start()**x).norm() << endl;
    }

    // setting error of macro iteration
    const double gradient = sigma->dot_product(*sigma) / sigma->size();

    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    if (gradient < thresh_) break;

    if (iter == max_iter_-1)
      throw runtime_error("Max iteration reached in the CASSCF macro interation.");
  }
  // ============================
  // macro iteration to here
  // ============================

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  fci_->update(coeff_);
  fci_->compute();
  fci_->compute_rdm12();
}


shared_ptr<const RotFile> CASBFGS::compute_denom(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {
  auto out = make_shared<RotFile>(nclosed_, nact_, nvirt_, false);
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


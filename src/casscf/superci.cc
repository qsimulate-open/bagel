//
// BAGEL - Parallel electron correlation program.
// Filename: superci.cc
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
#include <src/casscf/supercimicro.h>
#include <iostream>
#include <src/fci/fci.h>
#include <src/casscf/rotfile.h>
#include <src/math/davidson.h>
#include <src/molecule/hcore.h>
#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <src/math/bfgs.h>
#include <src/math/hpw_diis.h>

using namespace std;
using namespace bagel;

void SuperCI::compute() {
  // DIIS: will be turned on at iter = diis_start_ (>1),
  //       update log(U) where Cnow = Corig U. This is basically the same as the Hampel-Peterson-Werner
  //       paper on Brueckner CC
  shared_ptr<HPW_DIIS<Matrix>> diis;

  // BFGS: optional quasi-second-order MCSCF
//shared_ptr<BFGS<RotFile>> bfgs;

  // ============================
  // macro iteration from here
  // ============================
  double gradient = 1.0e100;
  mute_stdcout();
  Timer timer;
  for (int iter = 0; iter != max_iter_; ++iter) {

    if (iter >= diis_start_ && gradient < 1.0e-2 && diis == nullptr) {
      shared_ptr<Matrix> tmp = coeff_->copy();
      shared_ptr<Matrix> unit = make_shared<Matrix>(coeff_->mdim(), coeff_->mdim());
      unit->unit();
      diis = make_shared<HPW_DIIS<Matrix>>(10, tmp, unit);
    }

    // first perform CASCI to obtain RDMs
    if (iter) fci_->update(coeff_);
    Timer fci_time(0);
    fci_->compute();
    fci_->compute_rdm12();
    fci_time.tick_print("FCI and RDMs");
    // get energy
    energy_ = fci_->energy();

    // here make a natural orbitals and update the coefficients
    shared_ptr<Matrix> natorb = form_natural_orbs();

    auto grad = make_shared<RotFile>(nclosed_, nact_, nvirt_);

    // compute one-boedy operators
    shared_ptr<Matrix> f, fact, factp, gaa;
    shared_ptr<RotFile> denom;
    Timer onebody(0);
    one_body_operators(f, fact, factp, gaa, denom);
    onebody.tick_print("One body operators");

    // first, <proj|H|0> is computed
    grad->zero();
    // <a/i|H|0> = 2f_ai
    grad_vc(f, grad);
    // <a/r|H|0> = h_as d_sr + (as|tu)D_rs,tu = fact_ar
    grad_va(fact, grad);
    // <r/i|H|0> = 2f_ri - f^inact_is d_sr - 2(is|tu)P_rs,tu = 2f_ri - fact_ri
    grad_ca(f, fact, grad);

    // setting error of macro iteration
    gradient = grad->rms();
    if (gradient < thresh_) {
      rms_grad_ = gradient;
      break;
    }

    shared_ptr<const RotFile> cc;
    {
      Timer microiter_time(0);
      SuperCIMicro micro(shared_from_this(), grad, denom, f, fact, factp, gaa);
      micro.compute();
      cc = micro.cc();
      microiter_time.tick_print("Microiterations");
    }

    // unitary matrix
    shared_ptr<Matrix> rot = cc->unpack<Matrix>()->exp();
    // forcing rot to be unitary (usually not needed, though)
    rot->purify_unitary();

    if (diis == nullptr) {
      coeff_ = make_shared<const Coeff>(*coeff_ * *rot);
    } else {
      // including natorb.first to rot so that they can be processed at once
      shared_ptr<Matrix> tmp = rot->copy();
      dgemm_("N", "N", nact_, nbasis_, nact_, 1.0, natorb->data(), nact_, rot->element_ptr(nclosed_, 0), nbasis_, 0.0,
                                                                          tmp->element_ptr(nclosed_, 0), nbasis_);
      shared_ptr<const Matrix> tmp2 = tailor_rotation(tmp)->copy();
      shared_ptr<const Matrix> mcc = diis->extrapolate(tmp2);
      coeff_ = make_shared<const Coeff>(*mcc);
    }

    // synchronization
    mpi__->broadcast(const_pointer_cast<Coeff>(coeff_)->data(), coeff_->size(), 0);

    // print out...
    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      if (rms_grad_ > thresh_) cout << "    * The calculation did NOT converge. *    " << endl;
      cout << "    * Max iteration reached in the CASSCF macro interations. *     " << endl << endl;
    }
    mute_stdcout();

  }
  // ============================
  // macro iteration to here
  // ============================
  resume_stdcout();

  // block diagonalize coeff_ in nclosed and nvirt
  coeff_ = semi_canonical_orb();

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  fci_->update(coeff_);
  fci_->compute();
  fci_->compute_rdm12();
}


// rotate (within allowed rotations) the transformation matrix so that it is diagonal in each subblock
shared_ptr<Matrix> SuperCI::tailor_rotation(const shared_ptr<Matrix> seed) {

  shared_ptr<Matrix> out = seed->clone();
  for (int i = 0; i != nclosed_; ++i)
    for (int j = 0; j != nclosed_; ++j)
      out->element(j,i) = seed->element(j,i);
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      out->element(j+nclosed_,i+nclosed_) = seed->element(j+nclosed_,i+nclosed_);
  for (int i = 0; i != nvirt_; ++i)
    for (int j = 0; j != nvirt_; ++j)
      out->element(j+nocc_,i+nocc_) = seed->element(j+nocc_,i+nocc_);
  out->inverse();
  out->purify_unitary();
  *out = *seed * *out;

  return out;
}



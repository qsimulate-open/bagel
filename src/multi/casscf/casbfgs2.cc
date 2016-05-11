//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: casbfgs2.cc
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
#include <src/multi/casscf/qvec.h>
#include <src/scf/hf/fock.h>
#include <src/util/math/step_restrict_bfgs2.h>
#include <src/prop/hyperfine.h>

using namespace std;
using namespace bagel;

void CASBFGS2::compute() {

  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.

  shared_ptr<SRBFGS2<RotFile>> bfgs;
  shared_ptr<const Matrix> coeff_orig = coeff_->copy();
  auto arot = make_shared<RotFile>(nclosed_, nact_, nvirt_);

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;
  double oldenergy = 0.0;
  double oldgrad   = 1.0;

  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    if (nact_) {
      if (iter) fci_->update(coeff_);
      Timer fci_time(0);
      fci_->compute();
      fci_->compute_rdm12();
      fci_time.tick_print("FCI and RDMs");
      occup_ = VectorB(nact_);
      for (int i = 0; i != nact_; ++i)
        occup_[i] = fci_->rdm1_av()->element(i,i);
    }

    auto sigma = make_shared<RotFile>(nclosed_, nact_, nvirt_);

    // compute one-body operators
    Timer onebody;
    // * preparation
    const MatView ccoeff = coeff_->slice(0, nclosed_);
    // * core Fock operator
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<const Fock<1>>(geom_, hcore_, nullptr, ccoeff, /*store*/false, /*rhf*/true) : hcore_;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    // * active Fock operator
    // first make a weighted coefficient
    shared_ptr<Matrix> acoeff;
    if (nact_) {
      Matrix rdm1av(nact_, nact_);
      copy_n(fci_->rdm1_av()->data(), nact_*nact_, rdm1av.data());
      rdm1av.sqrt();
      rdm1av.scale(sqrt(1.0/2.0));
      acoeff = make_shared<Matrix>(coeff_->slice(nclosed_, nocc_) * rdm1av); 
    }
    // then make a AO density matrix
    shared_ptr<const Matrix> afock;
    if (nact_) {
      auto afockao = make_shared<Fock<1>>(geom_, hcore_, nullptr, acoeff, /*store*/false, /*rhf*/true);
      afock = make_shared<Matrix>(*coeff_ % (*afockao - *hcore_) * *coeff_);
    } else {
      afock = cfock->clone();
    }
    // * Q_xr = 2(xs|tu)P_rs,tu (x=general, mo)
    shared_ptr<const Qvec> qxr;
    if (nact_)
      qxr = make_shared<const Qvec>(coeff_->mdim(), nact_, coeff_, nclosed_, fci_, fci_->rdm2_av());

    // grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
    grad_vc(cfock, afock, sigma);
    // grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
    grad_va(cfock, qxr, sigma);
    // grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
    grad_ca(cfock, afock, qxr, sigma);

    // if this is the first time, set up the BFGS solver
    if (iter == 0) {
      shared_ptr<const RotFile> denom = compute_denom(cfock, afock, qxr);
      bfgs = make_shared<SRBFGS2<RotFile>>(denom);
    }
    onebody.tick_print("One body operators");

    // get energy
    if (nact_)
      energy_ = fci_->energy();
    else
      energy_ =  {(ccoeff % (*cfockao + *hcore_) * ccoeff).trace() + geom_->nuclear_repulsion()};

    // extrapolation using BFGS
    arot = bfgs->extrapolate(sigma, arot, blas::average(energy_));

    // restore the matrix from RotFile
    shared_ptr<const Matrix> amat = arot->unpack<Matrix>();
    shared_ptr<Matrix> expa = amat->exp(100);
    expa->purify_unitary();

    // updating coefficients
    coeff_ = make_shared<const Coeff>(*coeff_orig * *expa);
    // synchronization
    const_pointer_cast<Coeff>(coeff_)->synchronize();

    // setting error of macro iteration
    const double gradient = sigma->rms();
    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    if (gradient < thresh_) {
      rms_grad_ = gradient;
      cout << " " << endl;
      cout << "    * quasi-Newton optimization converged. *" << endl << endl;
      mute_stdcout();
      break;
    } else if (fabs(oldgrad - gradient) < thresh_ && fabs(blas::average(energy_)-oldenergy) < 1.0e-12) {
      rms_grad_ = gradient;
      only_energy_converged_ = true;
      cout << " " << endl;
      cout << "    * quasi-Newton optimization converged (only the average energy is converged) *" << endl << endl;
      mute_stdcout();
      break;
    }
    oldenergy = blas::average(energy_);
    oldgrad = gradient;

    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      if (rms_grad_ > thresh_) cout << "    * The calculation did NOT converge. *" << endl;
      cout << "    * Max iteration reached during the quasi-Newton optimization. *" << endl << endl;
    }
    mute_stdcout();
  }
  resume_stdcout();
  // ============================
  // macro iteration to here
  // ============================

  // block diagonalize coeff_ in nclosed and nvirt
  coeff_ = semi_canonical_orb();

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  if (nact_) {
    fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
  }

  // calculate the HFCCs
  if (do_hyperfine_ && !geom_->external() && nstate_ == 1) {
    HyperFine hfcc(geom_, spin_density(), fci_->det()->nspin(), "CASSCF");
    hfcc.compute();
  }
}

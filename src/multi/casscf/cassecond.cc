//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cassecond.cc
// Copyright (C) 2016 Toru Shiozaki
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

#include <src/scf/hf/fock.h>
#include <src/multi/casscf/qvec.h>
#include <src/multi/casscf/cassecond.h>

using namespace std;
using namespace bagel;

void CASSecond::compute() {
cout << "thresh " << thresh_ << endl;
cout << "max " << max_iter_ << endl;
  Timer timer;

  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    if (nact_) {
      if (iter) fci_->update(coeff_);
      Timer fci_time(0);
      fci_->compute();
      fci_->compute_rdm12();
      auto natorb = fci_->natorb_convert();
      coeff_ = update_coeff(coeff_, natorb.first);
      fci_time.tick_print("FCI and RDMs");
    }

    const MatView ccoeff = coeff_->slice(0, nclosed_);
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<const Fock<1>>(geom_, hcore_, nullptr, ccoeff, /*store*/false, /*rhf*/true) : hcore_;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    shared_ptr<const Matrix> afock;
    if (nact_) {
      shared_ptr<Matrix> acoeff = coeff_->slice_copy(nclosed_, nocc_);
      Matrix dkl(nact_, nact_);
      copy_n(fci_->rdm1_av()->data(), nact_*nact_, dkl.data());
      dkl.sqrt();
      dkl.scale(1.0/sqrt(2.0));
      auto afockao = make_shared<Fock<1>>(geom_, hcore_, nullptr, *acoeff * dkl, /*store*/false, /*rhf*/true);
      afock = make_shared<Matrix>(*coeff_ % (*afockao - *hcore_) * *coeff_);
    } else {
      afock = cfock->clone();
    }
    shared_ptr<const Qvec> qxr;
    if (nact_)
      qxr = make_shared<const Qvec>(coeff_->mdim(), nact_, coeff_, nclosed_, fci_, fci_->rdm2_av());

    auto sigma = make_shared<RotFile>(nclosed_, nact_, nvirt_);
    grad_vc(cfock, afock, sigma);
    grad_va(cfock, qxr, sigma);
    grad_ca(cfock, afock, qxr, sigma);

    sigma->print();

const double gradient = 0.0;
energy_ = {0.0, 0.0};
    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    if (gradient < thresh_) {
      cout << endl << "    * Second-order optimization converged. *   " << endl << endl;
      mute_stdcout();
      break;
    }

    if (iter == max_iter_-1)
      cout << endl << "    * Max iteration reached during the second optimization. *     " << endl << endl;
    mute_stdcout();
  }
  resume_stdcout();

  throw logic_error("reached the end of tmp implementaiton");
}


// grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
void CASSecond::grad_vc(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<RotFile> sigma) const {
  if (!nvirt_ || !nclosed_) return;
  double* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 4.0, cfock->element_ptr(nocc_,i), 1, target, 1);
    daxpy_(nvirt_, 4.0, afock->element_ptr(nocc_,i), 1, target, 1);
  }
}


// grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
// gamma is assumed to be diagonal
void CASSecond::grad_va(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> qxr, shared_ptr<RotFile> sigma) const {
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
void CASSecond::grad_ca(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<RotFile> sigma) const {
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


//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcassecond.cc
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

#include <src/scf/dhf/dfock.h>
#include <src/multi/zcasscf/zqvec.h>
#include <src/multi/zcasscf/zcassecond.h>

using namespace std;
using namespace bagel;

void ZCASSecond::compute() {
  assert(nvirt_ && nact_);
  Timer timer;

  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    {   
      if (iter) fci_->update(coeff_, /*restricted*/true);
      Timer fci_time(0);
      fci_->compute();
      fci_->compute_rdm12();
//    auto natorb = fci_->natorb_convert();
//    coeff_ = update_coeff(coeff_, natorb.first);
//    occup_ = natorb.second;
      fci_time.tick_print("FCI and RDMs");
      energy_ = fci_->energy();
    }

    shared_ptr<const ZMatrix> cfockao = fci_->jop()->core_fock();
    shared_ptr<const ZMatrix> afockao = compute_active_fock(coeff_->slice(nclosed_*2, nocc_*2), fci_->rdm1_av());
    shared_ptr<const ZMatrix> cfock = make_shared<ZMatrix>(*coeff_ % *cfockao * *coeff_);
    shared_ptr<const ZMatrix> afock = make_shared<ZMatrix>(*coeff_ % *afockao * *coeff_);
    shared_ptr<const ZMatrix> qxr = make_shared<ZQvec>(nbasis_*2, nact_, geom_, coeff_, coeff_->slice_copy(nclosed_*2,nocc_*2), nclosed_, fci_, gaunt_, breit_)->get_conjg();

    shared_ptr<const ZRotFile> grad = compute_gradient(cfock, afock, qxr);
throw logic_error("STOP!");

    // check gradient and break if converged
    const double gradient = grad->rms();
    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());
    mute_stdcout();
    if (gradient < thresh_) {
      resume_stdcout();
      cout << endl << "    * Second-order optimization converged. *   " << endl << endl;
      mute_stdcout();
      break;
    }
  }
}



shared_ptr<ZRotFile> ZCASSecond::compute_gradient(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr) const {
  auto sigma = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);
  shared_ptr<const ZMatrix> rdm1 = fci_->rdm1_av();
  if (nvirt_ && nclosed_) {
    complex<double>* target = sigma->ptr_vc();
    for (int i = 0; i != nclosed_*2; ++i, target += nvirt_*2) {
      zaxpy_(nvirt_*2, 1.0, cfock->element_ptr(nocc_*2, i), 1, target, 1);
      zaxpy_(nvirt_*2, 1.0, afock->element_ptr(nocc_*2, i), 1, target, 1);
    }
  }

  if (nvirt_ && nact_) {
    zgemm3m_("N", "T", nvirt_*2, nact_*2, nact_*2, 1.0, cfock->element_ptr(nocc_*2, nclosed_*2), cfock->ndim(), rdm1->data(), rdm1->ndim(), 0.0, sigma->ptr_va(), nvirt_*2);
    complex<double>* target = sigma->ptr_va();
    for (int i = 0; i != nact_*2; ++i, target += nvirt_*2) {
      zaxpy_(nvirt_*2, 1.0, qxr->element_ptr(nocc_*2, i), 1, target, 1);
    }
  }

  if (nclosed_ && nact_) {
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
  *sigma *= 2.0;
  return sigma;
}



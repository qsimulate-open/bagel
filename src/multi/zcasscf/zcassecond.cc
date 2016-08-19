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

#include <src/util/math/aughess.h>
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
      if (iter) fci_->update(coeff_);
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

    shared_ptr<ZRotFile> grad = compute_gradient(cfock, afock, qxr);

    // DEBUG CODE *************************
    zero_positronic_elements(grad); 
    shared_ptr<RelDFHalf> half, halfa;
    // DEBUG CODE *************************

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

    // compute denominator...
    shared_ptr<const ZRotFile> denom = compute_denom(cfock, afock, qxr, fci_->rdm1_av());

    AugHess<ZRotFile> solver(max_micro_iter_, grad);
    // initial trial vector
    shared_ptr<ZRotFile> trot = apply_denom(grad, denom, 0.001, 1.0);
    trot->normalize();

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      Timer mtimer;
      shared_ptr<const ZRotFile> sigma = compute_hess_trial(trot, half, halfa, cfock, afock, qxr);
      shared_ptr<const ZRotFile> residual;
      double lambda, epsilon, stepsize;
      tie(residual, lambda, epsilon, stepsize) = solver.compute_residual(trot, sigma);
      const double err = residual->norm() / lambda;
      resume_stdcout();
      cout << "         residual: " << setw(10) << setprecision(3) << scientific << err
           <<         " lambda  : " << setw(10) << setprecision(3) << scientific << lambda
           <<         " epsilon : " << setw(10) << setprecision(3) << scientific << epsilon
           <<         " stepsize: " << setw(10) << setprecision(3) << scientific << stepsize
           << setw(6) << fixed << setprecision(2) << mtimer.tick() << endl;
      mute_stdcout();
      if (err < max(thresh_micro_, stepsize*thresh_microstep_))
        break;

      trot = apply_denom(residual, denom, -epsilon, 1.0/lambda);
      for (int i = 0; i != 10; ++i) {
        const double norm = solver.orthog(trot);
        if (norm > 0.25) break;
      }
    }
grad->print();
throw logic_error("STOP!");
  }
}


shared_ptr<ZRotFile> ZCASSecond::apply_denom(shared_ptr<const ZRotFile> grad, shared_ptr<const ZRotFile> denom, const double shift, const double scale) const {
  shared_ptr<ZRotFile> out = grad->copy();
  for (int i = 0; i != out->size(); ++i)
    if (abs(denom->data(i)*scale+shift) > 1.0e-12)
      out->data(i) /= denom->data(i)*scale+shift;
  return out;
}


shared_ptr<ZRotFile> ZCASSecond::compute_hess_trial(shared_ptr<const ZRotFile> trot, shared_ptr<const RelDFHalf> half, shared_ptr<const RelDFHalf> halfa,
                                                    shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr) const {
  return nullptr;
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


// TODO this is an approximate denominator. We will replace this with exact diagonal.
shared_ptr<ZRotFile> ZCASSecond::compute_denom(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1) const {
  auto out = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);
  auto cfockd = make_shared<ZMatrix>(*cfock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2) * *rdm1);

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


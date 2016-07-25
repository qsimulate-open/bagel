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
      auto natorb = fci_->natorb_convert();
      coeff_ = update_coeff(coeff_, natorb.first);
      fci_time.tick_print("FCI and RDMs");
      energy_ = fci_->energy();
    }

    // coeff spaces and RDM1
    const MatView ccoeff = coeff_->slice(0, nclosed_);
    const MatView acoeff = coeff_->slice(nclosed_, nocc_);
    const MatView vcoeff = coeff_->slice(nocc_, nocc_+nvirt_);

    Matrix rdm1(nact_, nact_);
    copy_n(fci_->rdm1_av()->data(), nact_*nact_, rdm1.data());

    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<Fock<1>>(geom_, hcore_, nullptr, ccoeff, /*store*/true, /*rhf*/true) : hcore_;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    shared_ptr<const Matrix> afock;
    {
      Matrix dkl = rdm1;
      dkl.sqrt();
      dkl.scale(1.0/sqrt(2.0));
      auto afockao = make_shared<Fock<1>>(geom_, hcore_, nullptr, acoeff * dkl, /*store*/false, /*rhf*/true);
      afock = make_shared<Matrix>(*coeff_ % (*afockao - *hcore_) * *coeff_);
    }
    shared_ptr<const Qvec> qxr = make_shared<Qvec>(coeff_->mdim(), nact_, coeff_, nclosed_, fci_, fci_->rdm2_av());

    auto grad = make_shared<RotFile>(nclosed_, nact_, nvirt_);
    grad_vc(cfock, afock, grad);
    grad_va(cfock, qxr, grad);
    grad_ca(cfock, afock, qxr, grad);

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

    // half-transformed integrals (with JJ)
    shared_ptr<const DFHalfDist> half = nclosed_ ? dynamic_pointer_cast<const Fock<1>>(cfockao)->half()->apply_J() : nullptr;
    shared_ptr<const DFHalfDist> halfa = fci_->jop()->mo2e_1ext()->apply_JJ();

    // fock submatrices
    shared_ptr<const Matrix> fcaa = cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_);
    shared_ptr<const Matrix> faaa = afock->get_submatrix(nclosed_, nclosed_, nact_, nact_);
    shared_ptr<const Matrix> fcva = cfock->get_submatrix(nocc_, nclosed_, nvirt_, nact_);
    shared_ptr<const Matrix> fava = afock->get_submatrix(nocc_, nclosed_, nvirt_, nact_);
    shared_ptr<const Matrix> fcvv = cfock->get_submatrix(nocc_, nocc_, nvirt_, nvirt_);
    shared_ptr<const Matrix> favv = afock->get_submatrix(nocc_, nocc_, nvirt_, nvirt_);
    shared_ptr<const Matrix> fccc = nclosed_ ? cfock->get_submatrix(0, 0, nclosed_, nclosed_) : nullptr;
    shared_ptr<const Matrix> facc = nclosed_ ? afock->get_submatrix(0, 0, nclosed_, nclosed_) : nullptr;
    shared_ptr<const Matrix> fcca = nclosed_ ? cfock->get_submatrix(0, nclosed_, nclosed_, nact_) : nullptr;
    shared_ptr<const Matrix> faca = nclosed_ ? afock->get_submatrix(0, nclosed_, nclosed_, nact_) : nullptr;
    shared_ptr<const Matrix> fcvc = nclosed_ ? cfock->get_submatrix(nocc_, 0, nvirt_, nclosed_) : nullptr;
    shared_ptr<const Matrix> favc = nclosed_ ? afock->get_submatrix(nocc_, 0, nvirt_, nclosed_) : nullptr;

    // compute denominator...
    auto denom = grad->clone();
    {
      const Matrix fcd = *fcaa * rdm1;
      for (int i = 0; i != nact_; ++i)
        for (int j = 0; j != nclosed_; ++j)
          denom->ele_ca(j, i) += 4.0 * ((*fcaa)(i, i) + (*faaa)(i, i)) - 4.0 * ((*fccc)(j, j) + (*facc)(j, j)) - 2.0 * fcd(i, i) + 2.0 * (*fccc)(j, j) * rdm1(i, i);
      for (int i = 0; i != nclosed_; ++i)
        for (int j = 0; j != nvirt_; ++j)
          denom->ele_vc(j, i) += 4.0 * ((*fcvv)(j, j) + (*favv)(j, j)) - 4.0 * ((*fccc)(i, i) + (*facc)(i, i));
      for (int i = 0; i != nact_; ++i)
        for (int j = 0; j != nvirt_; ++j)
          denom->ele_va(j, i) += 2.0 * (*fcvv)(j, j) * rdm1(i, i) - 2.0 * fcd(i, i);
    }
    denom->print();

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      Timer mtimer;
      // trial vector
      auto trot = grad->clone();
trot->ax_plus_y_vc(1.0, *grad->vc_mat());
trot->ax_plus_y_va(1.0, *grad->va_mat());
trot->ax_plus_y_ca(1.0, *grad->ca_mat());
trot->print();
      shared_ptr<const Matrix> va = trot->va_mat();
      shared_ptr<const Matrix> ca = nclosed_ ? trot->ca_mat() : nullptr;
      shared_ptr<const Matrix> vc = nclosed_ ? trot->vc_mat() : nullptr;
      // sigma vector
      auto sigma = grad->clone();

      // lambda for computing g(D)
      auto compute_gd = [&,this](shared_ptr<const DFHalfDist> halft, shared_ptr<const DFHalfDist> halfjj, const MatView pcoeff) {
        shared_ptr<const Matrix> pcoefft = make_shared<Matrix>(pcoeff)->transpose();
        shared_ptr<Matrix> gd = geom_->df()->compute_Jop(halft, pcoefft);
        shared_ptr<Matrix> ex0 = halfjj->form_2index(halft, 1.0);
        ex0->symmetrize();
        gd->ax_plus_y(-0.5, ex0);
        return gd;
      };

      // g(t_vc) operator and g(t_ac) operator
      if (nclosed_) {
        const Matrix tcoeff = vcoeff * *vc + acoeff * *ca->transpose();
        auto halft = geom_->df()->compute_half_transform(tcoeff);
        const Matrix gt = *compute_gd(halft, half, ccoeff);
        sigma->ax_plus_y_ca(32.0, ccoeff % gt * acoeff);
        sigma->ax_plus_y_vc(32.0, vcoeff % gt * ccoeff);
        sigma->ax_plus_y_va(16.0, vcoeff % gt * acoeff * rdm1);
        sigma->ax_plus_y_ca(-16.0, ccoeff % gt * acoeff * rdm1);
      }
      // g(t_va - t_ca)
      const Matrix tcoeff = nclosed_ ? (vcoeff * *va - ccoeff * *ca) : vcoeff * *va;
      shared_ptr<const DFHalfDist> halfta = geom_->df()->compute_half_transform(tcoeff);
      if (nclosed_) {
        shared_ptr<DFHalfDist> halftad = halfta->copy();
        halftad->rotate_occ(make_shared<Matrix>(rdm1));
        const Matrix gt = *compute_gd(halftad, halfa, acoeff);
        sigma->ax_plus_y_ca(16.0, ccoeff % gt * acoeff);
        sigma->ax_plus_y_vc(16.0, vcoeff % gt * ccoeff);
      }
      // terms with Qvec
      {
        shared_ptr<const Matrix> qaa = qxr->cut(nclosed_, nocc_);
        shared_ptr<const Matrix> qva = qxr->cut(nocc_, nocc_+nvirt_);
        sigma->ax_plus_y_va(-2.0, *va ^ *qaa);
        sigma->ax_plus_y_va(-2.0, *va * *qaa);
        if (nclosed_) {
          shared_ptr<const Matrix> qca = qxr->cut(0, nclosed_);
          sigma->ax_plus_y_vc(-2.0, *va ^ *qca);
          sigma->ax_plus_y_va(-2.0, *vc * *qca);
          sigma->ax_plus_y_ca(-2.0, *vc % *qva);
          sigma->ax_plus_y_vc(-2.0, *qva ^ *ca);
          sigma->ax_plus_y_ca(-2.0, *ca ^ *qaa);
          sigma->ax_plus_y_ca(-2.0, *ca * *qaa);
        }
      }
      // compute Q' and Q''
      {
        shared_ptr<const DFFullDist> fullaa = halfa->compute_second_transform(acoeff);
        shared_ptr<DFFullDist> fullta = halfta->compute_second_transform(acoeff);
        shared_ptr<const DFFullDist> fulltas = fullta->swap();
        fullta->ax_plus_y(1.0, fulltas);
        shared_ptr<const DFFullDist> fullaaD = fullaa->apply_2rdm(*fci_->rdm2_av());
        shared_ptr<const DFFullDist> fulltaD = fullta->apply_2rdm(*fci_->rdm2_av());
        shared_ptr<const Matrix> qp  = halfa->form_2index(fulltaD, 1.0);
        shared_ptr<const Matrix> qpp = halfta->form_2index(fullaaD, 1.0);

        sigma->ax_plus_y_va( 4.0, vcoeff % (*qp + *qpp));
        if (nclosed_)
          sigma->ax_plus_y_ca(-4.0, ccoeff % (*qp + *qpp));
      }

      // next 1-electron contribution...
      {
        sigma->ax_plus_y_va( 4.0, *fcvv * *va * rdm1);
        sigma->ax_plus_y_va(-2.0, *va * (rdm1 * *fcaa + *fcaa * rdm1));
        if (nclosed_) {
          sigma->ax_plus_y_ca( 8.0, *ca * (*fcaa + *faaa));
          sigma->ax_plus_y_ca( 8.0, *vc % (*fcva + *fava));
          sigma->ax_plus_y_vc(-8.0, *vc * (*fccc + *facc));
          sigma->ax_plus_y_va(-4.0, *vc * (*fcca + *faca));
          sigma->ax_plus_y_vc(-4.0, *va ^ (*fcca + *faca));
          sigma->ax_plus_y_ca(-2.0, *ca * (rdm1 * *fcaa + *fcaa * rdm1));
          sigma->ax_plus_y_vc( 8.0, (*fcvv + *favv) * *vc);
          sigma->ax_plus_y_ca(-8.0, (*fccc + *facc) * *ca);
          sigma->ax_plus_y_va( 4.0, (*fcvc + *favc) * *ca);
          sigma->ax_plus_y_ca( 4.0, (*fcvc + *favc) % *va);
          sigma->ax_plus_y_vc( 4.0, (*fcva + *fava) ^ *ca);
          sigma->ax_plus_y_vc( 4.0, (*fcva + *fava) ^ *ca);
          sigma->ax_plus_y_ca( 4.0, *fccc * *ca * rdm1);
          sigma->ax_plus_y_ca(-4.0, *fcvc % *va * rdm1);
          sigma->ax_plus_y_va(-4.0, *fcvc * *ca * rdm1);
          sigma->ax_plus_y_vc(-2.0, *fcva * rdm1 ^ *ca);
          sigma->ax_plus_y_vc(-2.0, *va * rdm1 ^ *fcca);
          sigma->ax_plus_y_ca(-2.0, *vc % *fcva * rdm1);
          sigma->ax_plus_y_va(-2.0, *vc * *fcca * rdm1);
        }
      }

      sigma->scale(0.5);

sigma->print();
      cout << "         " << mtimer.tick() << endl;
break;
    }

    resume_stdcout();
    if (iter == max_iter_-1)
      cout << endl << "    * Max iteration reached during the second optimization. *     " << endl << endl;
    mute_stdcout();
break;
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


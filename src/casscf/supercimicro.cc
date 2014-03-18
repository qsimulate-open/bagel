//
// BAGEL - Parallel electron correlation program.
// Filename: superci_micro.cc
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

#include <src/math/pairfile.h>
#include <src/math/bfgs.h>
#include <src/casscf/supercimicro.h>

using namespace std;
using namespace bagel;


void SuperCIMicro::compute() {

  using SCIData = PairFile<Matrix,RotFile>;

  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
  const int nvirt = casscf_->nvirt();
  // a pair of 1*1 matrix and RotFile
  DavidsonDiag<SCIData> davidson(1, casscf_->max_micro_iter());

  // current coefficient
  auto cc0    = make_shared<Matrix>(1,1,true);
  auto cc1    = make_shared<RotFile>(nclosed, nact, nvirt);
  // BFGS initialization
  auto mbfgs = make_shared<BFGS<SCIData>>(make_shared<SCIData>(cc0, make_shared<RotFile>(*denom_)));

  for (int miter = 0; miter != casscf_->max_micro_iter(); ++miter) {
    Timer mtimer;

    shared_ptr<RotFile> sigma1;
    auto sigma0 = cc0->clone();

    if (miter != 0) {
      sigma1 = form_sigma(cc1);
      // projection to reference
      (*cc0)   (0,0) = 0.0;
      (*sigma0)(0,0) = grad_->dot_product(*cc1);
    } else {
      sigma1 = grad_->copy();
      (*cc0)   (0,0) = 1.0;
      (*sigma0)(0,0) = 0.0;
    }

    cc1->synchronize();
    sigma1->synchronize();

    // enters davidson iteration
    auto ccp    = make_shared<SCIData>(cc0->copy(), cc1->copy());
    auto sigmap = make_shared<SCIData>(sigma0->copy(), sigma1->copy());
    const double mic_energy = davidson.compute(ccp, sigmap);

    // residual vector and error
    shared_ptr<SCIData> residual = davidson.residual().front();
    const double error = residual->rms();

    if (miter == 0) cout << endl << "     == micro iteration == " << endl;
    cout << setw(10) << miter << "   " << setw(20) << setprecision(12) << mic_energy << " "
         << setw(10) << scientific << setprecision(2) << error << fixed << " " << mtimer.tick() << endl;

    if (error < casscf_->thresh_micro()) { cout << endl; break; }
    if (miter+1 == casscf_->max_micro_iter()) throw runtime_error("max_micro_iter_ is reached in CASSCF");

    // update cc0 and cc1
    cc1 = mbfgs->extrapolate(residual, davidson.civec().front())->second();
    cc1->normalize();
    cc0 = cc0->clone();
  }

  // rotation parameters
  shared_ptr<const SCIData> result = davidson.civec().front();
  const double cref = result->first()->element(0,0);
  shared_ptr<RotFile> tmp = result->second()->copy();
  *tmp *= 1.0/cref;
  cc_ = tmp;
}


std::shared_ptr<RotFile> SuperCIMicro::form_sigma(std::shared_ptr<const RotFile> cc) const {

  auto sigma = cc->clone();
  // equation 21d
  sigma_at_at_(cc, sigma);
  // equation 21e
  sigma_ai_ai_(cc, sigma);
  // equation 21f // note a typo!
  sigma_at_ai_(cc, sigma);
  // equation 21b
  sigma_ai_ti_(cc, sigma);
  // equation 21a
  sigma_ti_ti_(cc, sigma);

  return sigma;
}


// sigma_at_at = delta_ab Gtu/sqrt(nt nu) + delta_tu Fab
void SuperCIMicro::sigma_at_at_(shared_ptr<const RotFile> cc, shared_ptr<RotFile> sigma) const {
  const int nact = casscf_->nact();
  const int nvirt = casscf_->nvirt();
  const int nocc = casscf_->nocc();
  const int nbasis = casscf_->nbasis();
  if (!nact || !nvirt) return;

  shared_ptr<Matrix> gtup = gaa_->copy();
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nact; ++j) {
#if 0
      const double fac = (occup_[i]*occup_[j] > occup_thresh) ? 1.0/std::sqrt(occup_[i]*occup_[j]) : 0.0;
#else
      const double fac = 1.0;
#endif
      gtup->element(j,i) *= fac;
    }
  }
  dgemm_("N", "N", nvirt, nact, nact, 1.0, cc->ptr_va(), nvirt, gtup->data(), nact, 1.0, sigma->ptr_va(), nvirt);
  dgemm_("N", "N", nvirt, nact, nvirt, 1.0, fock_->element_ptr(nocc, nocc), nbasis, cc->ptr_va(), nvirt, 1.0, sigma->ptr_va(), nvirt);
}


// sigma_ai_ai = delta_ij F_ab - delta_ab F_ij
void SuperCIMicro::sigma_ai_ai_(shared_ptr<const RotFile> cc, shared_ptr<RotFile> sigma) const {
  const int nclosed = casscf_->nclosed();
  const int nvirt = casscf_->nvirt();
  const int nocc = casscf_->nocc();
  const int nbasis = casscf_->nbasis();
  if (!nclosed || !nvirt) return;

  dgemm_("N", "N", nvirt, nclosed, nclosed, -1.0, cc->ptr_vc(), nvirt, fock_->data(), nbasis, 1.0, sigma->ptr_vc(), nvirt);
  dgemm_("N", "N", nvirt, nclosed, nvirt, 1.0, fock_->element_ptr(nocc, nocc), nbasis, cc->ptr_vc(), nvirt, 1.0, sigma->ptr_vc(), nvirt);
}


// sigma_at_ai = -delta_ab Fact_ti sqrt(nt/2)
void SuperCIMicro::sigma_at_ai_(shared_ptr<const RotFile> cc, shared_ptr<RotFile> sigma) const {
  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
  const int nvirt = casscf_->nvirt();
  if (!nact || !nvirt || !nclosed) return;

  Matrix tmp(nclosed, nact);
  tmp.zero();
  for (int i = 0; i != nact; ++i) {
    const double fac = -std::sqrt(0.5*casscf_->occup(i));
    daxpy_(nclosed, fac, fockact_->element_ptr(0,i), 1, tmp.element_ptr(0,i), 1);
  }
  dgemm_("N", "N", nvirt, nact, nclosed, 1.0, cc->ptr_vc(), nvirt, tmp.data(), nclosed, 1.0, sigma->ptr_va(), nvirt);
  dgemm_("N", "T", nvirt, nclosed, nact, 1.0, cc->ptr_va(), nvirt, tmp.data(), nclosed, 1.0, sigma->ptr_vc(), nvirt);
}


// sigma_ai_ti = sqrt((2-nt)/2) Fact_at
void SuperCIMicro::sigma_ai_ti_(shared_ptr<const RotFile> cc, shared_ptr<RotFile> sigma) const {
  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
  const int nvirt = casscf_->nvirt();
  const int nocc = casscf_->nocc();
  if (!nact || !nvirt || !nclosed) return;

  Matrix tmp(nvirt, nact);
  tmp.zero();
  for (int i = 0; i != nact; ++i) {
    const double fac = std::sqrt(1.0-0.5*casscf_->occup(i));
    daxpy_(nvirt, fac, fockact_->element_ptr(nocc,i), 1, tmp.element_ptr(0,i), 1);
  }
  dgemm_("T", "N", nclosed, nact, nvirt, 1.0, cc->ptr_vc(), nvirt, tmp.data(), nvirt, 1.0, sigma->ptr_ca(), nclosed);
  dgemm_("N", "T", nvirt, nclosed, nact, 1.0, tmp.data(), nvirt, cc->ptr_ca(), nclosed, 1.0, sigma->ptr_vc(), nvirt);

}


// sigma_ti_ti = - delta_ij ((2-nt-nu)Fact_tu - G_tu)/sqrt((2-nt)(2-nu)) - delta_tu f_ij
void SuperCIMicro::sigma_ti_ti_(shared_ptr<const RotFile> cc, shared_ptr<RotFile> sigma) const {
  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
  const int nbasis = casscf_->nbasis();
  if (!nact || !nclosed) return;
  Matrix tmp(nact, nact);
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nact; ++j) {
      const double fac = ((2.0-casscf_->occup(i))*(2.0-casscf_->occup(j)) > occup_thresh) ? 1.0/std::sqrt((2.0-casscf_->occup(i))*(2.0-casscf_->occup(j))) : 0.0;
      tmp(j,i) = -((2.0 - casscf_->occup(j) - casscf_->occup(i)) * fockactp_->element(j,i) - gaa_->element(j,i)) * fac;
    }
  }
  dgemm_("N", "N", nclosed, nact, nact, 1.0, cc->ptr_ca(), nclosed, tmp.data(), nact, 1.0, sigma->ptr_ca(), nclosed);
  dgemm_("N", "N", nclosed, nact, nclosed, -1.0, fock_->data(), nbasis, cc->ptr_ca(), nclosed, 1.0, sigma->ptr_ca(), nclosed);
}


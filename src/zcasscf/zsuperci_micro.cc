//
// BAGEL - Parallel electron correlation program.
// Filename: zsuperci_micro.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson Bates <jefferson.bates@northwestern.edu>
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
#include <src/zcasscf/zsupercimicro.h>
#include <src/zcasscf/zsuperci.h>

using namespace std;
using namespace bagel;

void ZSuperCIMicro::compute() {

  using ZSCIData = PairFile<ZMatrix,ZRotFile>;

  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
#ifdef BOTHSPACES
  const int nvirt = casscf_->nvirt();
#else
  const int nvirt = casscf_->nvirtnr();
#endif
  DavidsonDiag<ZSCIData, ZMatrix> davidson(1, casscf_->max_micro_iter());

  // current coefficient
  auto cc0    = make_shared<ZMatrix>(1,1,true);
  auto cc1    = make_shared<ZRotFile>(nclosed*2, nact*2, nvirt*2);
  // BFGS initialization
  auto mbfgs = make_shared<BFGS<ZSCIData>>(make_shared<ZSCIData>(cc0, make_shared<ZRotFile>(*denom_)));

  for (int miter = 0; miter != casscf_->max_micro_iter(); ++miter) {
    Timer mtimer;

    shared_ptr<ZRotFile> sigma1;
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

    // enters davidson iteration
    auto ccp    = make_shared<ZSCIData>(cc0->copy(), cc1->copy());
    auto sigmap = make_shared<ZSCIData>(sigma0->copy(), sigma1->copy());
//    ccp->synchronize();
//    sigmap->synchronize();
    const double mic_energy = davidson.compute(ccp, sigmap);

    // residual vector and error
    shared_ptr<ZSCIData> residual = davidson.residual().front();
//    residual->synchronize();
    const double error = residual->rms();
    assert(isnormal(error)); // check for nan's

    if (miter == 0) cout << endl << "     == micro iteration == " << endl;
    cout << setw(10) << miter << "   " << setw(20) << setprecision(12) << mic_energy << " "
         << setw(10) << scientific << setprecision(2) << error << fixed << " " << mtimer.tick() << endl;

    if (error < casscf_->thresh_micro()) { cout << endl; break; }
//    if (miter+1 == casscf_->max_micro_iter()) throw runtime_error("max_micro_iter_ is reached in CASSCF");

    // update cc0 and cc1
    cc1 = mbfgs->extrapolate(residual, davidson.civec().front())->second();
    cc1->normalize();
    cc0 = cc0->clone();
  }

  // rotation parameters
  shared_ptr<const ZSCIData> result = davidson.civec().front();
  const complex<double> cref = result->first()->element(0,0);
  shared_ptr<ZRotFile> tmp = result->second()->copy();
  *tmp *= 1.0/cref;
//  tmp->synchronize();
  cc_ = tmp;
}


std::shared_ptr<ZRotFile> ZSuperCIMicro::form_sigma(std::shared_ptr<const ZRotFile> cc) const {

  auto sigma = cc->clone();
  // equation 21f // note a typo!
  sigma_at_at_(cc, sigma);
  // equation 21d
  sigma_ai_ai_(cc, sigma);
  // equation 21e
  sigma_at_ai_(cc, sigma);
  // equation 21b
  sigma_ai_ti_(cc, sigma);
  // equation 21a
  sigma_ti_ti_(cc, sigma);

  return sigma;
}


// sigma_at_at = 2 * delta_ab Gtu/sqrt(nt nu) + delta_tu Fab  : TODO factor 2 needed here to reproduce non-rel limit
// TODO : check why normalization factor is commented out
void ZSuperCIMicro::sigma_at_at_(shared_ptr<const ZRotFile> cc, shared_ptr<ZRotFile> sigma) const {
  const int nact = casscf_->nact();
#ifdef BOTHSPACES
  const int nvirt = casscf_->nvirt();
  const int nbasis = casscf_->nbasis();
#else
  const int nvirt = casscf_->nvirtnr();
  const int nbasis = casscf_->nbasis()/2;
#endif
  const int nocc = casscf_->nocc();
  if (!nact || !nvirt) return;

  shared_ptr<ZMatrix> gtup = gaa_->copy();
  for (int i = 0; i != nact*2; ++i) {
    for (int j = 0; j != nact*2; ++j) {
#if 0
      const double fac = (casscf_->occup(i)*casscf_->occup(j) > zoccup_thresh) ? 1.0/std::sqrt(casscf_->occup(i)*casscf_->occup(j)) : 0.0;
#else
      const double fac = 2.0;
#endif
      gtup->element(j,i) *= fac;
    }
  }
  zgemm3m_("N", "N", nvirt*2, nact*2, nact*2, 1.0, cc->ptr_va(), nvirt*2, gtup->data(), nact*2, 1.0, sigma->ptr_va(), nvirt*2);
  zgemm3m_("N", "N", nvirt*2, nact*2, nvirt*2, 1.0, fock_->element_ptr(nocc*2, nocc*2), nbasis*2, cc->ptr_va(), nvirt*2, 1.0, sigma->ptr_va(), nvirt*2);
}


// sigma_ai_ai = delta_ij F_ab - delta_ab F_ij
void ZSuperCIMicro::sigma_ai_ai_(shared_ptr<const ZRotFile> cc, shared_ptr<ZRotFile> sigma) const {
  const int nclosed = casscf_->nclosed();
#ifdef BOTHSPACES
  const int nvirt = casscf_->nvirt();
  const int nbasis = casscf_->nbasis();
#else
  const int nvirt = casscf_->nvirtnr();
  const int nbasis = casscf_->nbasis()/2;
#endif
  const int nocc = casscf_->nocc();
  if (!nclosed || !nvirt) return;

  zgemm3m_("N", "N", nvirt*2, nclosed*2, nclosed*2, -1.0, cc->ptr_vc(), nvirt*2, fock_->data(), nbasis*2, 1.0, sigma->ptr_vc(), nvirt*2);
  zgemm3m_("N", "N", nvirt*2, nclosed*2, nvirt*2, 1.0, fock_->element_ptr(nocc*2, nocc*2), nbasis*2, cc->ptr_vc(), nvirt*2, 1.0, sigma->ptr_vc(), nvirt*2);
}


// sigma_at_ai = -delta_ab Fact_ti sqrt(nt)*2 TODO : not including sqrt(1/2) ; factor 2 needed to reproduce non-rel limit
void ZSuperCIMicro::sigma_at_ai_(shared_ptr<const ZRotFile> cc, shared_ptr<ZRotFile> sigma) const {
  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
#ifdef BOTHSPACES
  const int nvirt = casscf_->nvirt();
#else
  const int nvirt = casscf_->nvirtnr();
#endif
  if (!nact || !nvirt || !nclosed) return;

  ZMatrix tmp(nclosed*2, nact*2);
  tmp.zero();
  for (int i = 0; i != nact*2; ++i) {
    const double fac = -2.0 * std::sqrt(casscf_->occup(i)); // TODO check on a factor of 0.5
    zaxpy_(nclosed*2, fac, fockact_->element_ptr(0,i), 1, tmp.element_ptr(0,i), 1);
  }
  zgemm3m_("N", "N", nvirt*2, nact*2, nclosed*2, 1.0, cc->ptr_vc(), nvirt*2, tmp.data(), nclosed*2, 1.0, sigma->ptr_va(), nvirt*2);
  zgemm3m_("N", "C", nvirt*2, nclosed*2, nact*2, 1.0, cc->ptr_va(), nvirt*2, tmp.data(), nclosed*2, 1.0, sigma->ptr_vc(), nvirt*2);
}


// sigma_ai_ti = sqrt((1-nt))*2* Fact_at // TODO check normalization factors ; factor 2 needed to recover non-rel limit
void ZSuperCIMicro::sigma_ai_ti_(shared_ptr<const ZRotFile> cc, shared_ptr<ZRotFile> sigma) const {
  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
#ifdef BOTHSPACES
  const int nvirt = casscf_->nvirt();
#else
  const int nvirt = casscf_->nvirtnr();
#endif
  const int nocc = casscf_->nocc();
  if (!nact || !nvirt || !nclosed) return;

  ZMatrix tmp(nvirt*2, nact*2);
  tmp.zero();
  for (int i = 0; i != nact*2; ++i) {
    const double fac = 1.0 - casscf_->occup(i) > zoccup_thresh ? std::sqrt(1.0-casscf_->occup(i)) : 0.0;
    zaxpy_(nvirt*2, fac*2.0, fockact_->element_ptr(nocc*2,i), 1, tmp.element_ptr(0,i), 1);
  }
  zgemm3m_("C", "N", nclosed*2, nact*2, nvirt*2, 1.0, cc->ptr_vc(), nvirt*2, tmp.data(), nvirt*2, 1.0, sigma->ptr_ca(), nclosed*2); // TODO : check "C" vs "T" here
  zgemm3m_("N", "C", nvirt*2, nclosed*2, nact*2, 1.0, tmp.data(), nvirt*2, cc->ptr_ca(), nclosed*2, 1.0, sigma->ptr_vc(), nvirt*2);
}


// sigma_ti_ti = - delta_ij ((2-nt-nu)Fact_tu - G_tu)/sqrt((2-nt)(2-nu)) - delta_tu f_ij // TODO : check normalization factors
void ZSuperCIMicro::sigma_ti_ti_(shared_ptr<const ZRotFile> cc, shared_ptr<ZRotFile> sigma) const {
  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
#ifdef BOTHSPACES
  const int nbasis = casscf_->nbasis();
#else
  const int nbasis = casscf_->nbasis()/2;
#endif
  if (!nact || !nclosed) return;
  ZMatrix tmp(nact*2, nact*2);
  for (int i = 0; i != nact*2; ++i) {
    for (int j = 0; j != nact*2; ++j) {
      const double fac = ((1.0-casscf_->occup(i))*(1.0-casscf_->occup(j)) > zoccup_thresh) ? 1.0/std::sqrt((1.0-casscf_->occup(i))*(1.0-casscf_->occup(j))) : 0.0;
      tmp(j,i) = -((1.0 - casscf_->occup(j) - casscf_->occup(i)) * fockactp_->element(j,i) - gaa_->element(j,i)) * fac;
    }
  }
  zgemm3m_("N", "N", nclosed*2, nact*2, nact*2, 1.0, cc->ptr_ca(), nclosed*2, tmp.data(), nact*2, 1.0, sigma->ptr_ca(), nclosed*2);
  zgemm3m_("N", "N", nclosed*2, nact*2, nclosed*2, -1.0, fock_->data(), nbasis*2, cc->ptr_ca(), nclosed*2, 1.0, sigma->ptr_ca(), nclosed*2);
}

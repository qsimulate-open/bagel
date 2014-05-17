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

using namespace std;
using namespace bagel;

void ZSuperCIMicro::compute() {

  using ZSCIData = PairFile<ZMatrix,ZRotFile>;

  const int nclosed = casscf_->nclosed();
  const int nact = casscf_->nact();
  const int nvirt = casscf_->nvirt();
  // a pair of 1*1 matrix and RotFile
  DavidsonDiag<ZSCIData> davidson(1, casscf_->max_micro_iter());

  // current coefficient
  auto cc0    = make_shared<ZMatrix>(1,1,true);
  auto cc1    = make_shared<ZRotFile>(nclosed, nact, nvirt);
  // BFGS initialization
  auto mbfgs = make_shared<BFGS<ZSCIData>>(make_shared<ZSCIData>(cc0, make_shared<ZRotFile>(*denom_)));

  for (int miter = 0; miter != casscf_->max_micro_iter(); ++miter) {
    Timer mtimer;

    shared_ptr<ZRotFile> sigma1;
    auto sigma0 = cc0->clone();

//    if (miter != 0) {
//      sigma1 = form_sigma(cc1);
//      // projection to reference
//      (*cc0)   (0,0) = 0.0;
//      (*sigma0)(0,0) = grad_->dot_product(*cc1);
//    } else {
//      sigma1 = grad_->copy();
//      (*cc0)   (0,0) = 1.0;
//      (*sigma0)(0,0) = 0.0;
//    }
//
//    // enters davidson iteration
//    auto ccp    = make_shared<SCIData>(cc0->copy(), cc1->copy());
//    auto sigmap = make_shared<SCIData>(sigma0->copy(), sigma1->copy());
//    ccp->synchronize();
//    sigmap->synchronize();
//    const double mic_energy = davidson.compute(ccp, sigmap);
//
//    // residual vector and error
//    shared_ptr<SCIData> residual = davidson.residual().front();
//    residual->synchronize();
//    const double error = residual->rms();
//
//    if (miter == 0) cout << endl << "     == micro iteration == " << endl;
//    cout << setw(10) << miter << "   " << setw(20) << setprecision(12) << mic_energy << " "
//         << setw(10) << scientific << setprecision(2) << error << fixed << " " << mtimer.tick() << endl;
//
//    if (error < casscf_->thresh_micro()) { cout << endl; break; }
//    if (miter+1 == casscf_->max_micro_iter()) throw runtime_error("max_micro_iter_ is reached in CASSCF");
//
//    // update cc0 and cc1
//    cc1 = mbfgs->extrapolate(residual, davidson.civec().front())->second();
//    cc1->normalize();
//    cc0 = cc0->clone();
  }
//
//  // rotation parameters
//  shared_ptr<const SCIData> result = davidson.civec().front();
//  const double cref = result->first()->element(0,0);
//  shared_ptr<RotFile> tmp = result->second()->copy();
//  *tmp *= 1.0/cref;
//  tmp->synchronize();
//  cc_ = tmp;
}




// sigma_at_at = delta_ab Gtu/sqrt(nt nu) + delta_tu Fab 
// TODO : check why normalization factor is commented out
void ZSuperCIMicro::sigma_at_at_(shared_ptr<const ZRotFile> cc, shared_ptr<ZRotFile> sigma) const {
  const int nact = casscf_->nact();
  const int nvirt = casscf_->nvirt();
  const int nocc = casscf_->nocc();
  const int nbasis = casscf_->nbasis();
  if (!nact || !nvirt) return;

  shared_ptr<ZMatrix> gtup = gaa_->copy();
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
  zgemm3m_("N", "N", nvirt, nact, nact, 1.0, cc->ptr_va(), nvirt, gtup->data(), nact, 1.0, sigma->ptr_va(), nvirt);
  zgemm3m_("N", "N", nvirt, nact, nvirt, 1.0, fock_->element_ptr(nocc, nocc), nbasis, cc->ptr_va(), nvirt, 1.0, sigma->ptr_va(), nvirt);
}

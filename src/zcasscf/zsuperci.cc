//
// BAGEL - Parallel electron correlation program.
// Filename: zuperci.cc
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


#include <src/zcasscf/zqvec.h>
#include <src/math/step_restrict_bfgs.h>
#include <src/rel/dfock.h>
#include <src/zcasscf/zsuperci.h>
#include <src/zcasscf/zsupercimicro.h>
#include <src/rel/reloverlap.h>

using namespace std;
using namespace bagel;

void ZSuperCI::compute() {
    // TODO : DIIS needs to be introduced eventually

  // ============================
  // macro iteration from here
  // ============================
   Timer timer;

  // intialize coefficients
  init_kramers_coeff();

  cout << setprecision(8) << " kramers restricted re part rms = " << coeff_->get_real_part()->rms() << endl;
  cout << setprecision(8) << " kramers restricted im part rms = " << coeff_->get_imag_part()->rms() << endl;

  if (nact_)
    fci_->update(coeff_);

  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    if (nact_) {
      mute_stdcout(/*fci*/true);
      if (iter) fci_->update(coeff_);
      cout << " Executing FCI calculation in Cycle " << iter << endl;
      fci_->compute();
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      energy_.push_back((fci_->energy())[0]);
      resume_stdcout();
    }

    auto grad = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);

    // compute one-boedy operators
    shared_ptr<ZMatrix> f, fact, factp, gaa;
    shared_ptr<ZRotFile> denom;
    one_body_operators(f, fact, factp, gaa, denom);

    // first, <proj|H|0> is computed
    grad->zero();
    // <a/i|H|0> = f_ai
    grad_vc(f, grad);
    // <a/r|H|0> = h_as d_sr + (as|tu)D_rs,tu = fact_ar
    grad_va(fact, grad);
    // <r/i|H|0> = f_ri - f^inact_is d_sr - (is|tu)P_rs,tu = f_ri - fact_ri
    grad_ca(f, fact, grad);

     // setting error of macro iteration
     auto gradient = grad->rms();
     if (gradient < thresh_) break;

  // ============================
     // Micro-iterations go here
  // ============================
    shared_ptr<const ZRotFile> cc;
    {
      ZSuperCIMicro micro(shared_from_this(), grad, denom, f, fact, factp, gaa);
      micro.compute();
      cc = micro.cc();
    }

   // orbital rotation matrix
    shared_ptr<ZMatrix> amat = cc->unpack<ZMatrix>();
    kramers_adapt(amat, nvirt_);

   // multiply multiply -i to make amat hermite (will be compensated), then make Exp(Kappa)
   *amat *= 1.0 * complex<double>(0.0, -1.0);
   unique_ptr<double[]> teig(new double[amat->ndim()]);
   amat->diagonalize(teig.get());
   auto amat_sav = amat->copy();
   for (int i = 0; i != amat->ndim(); ++i) {
     complex<double> ex = exp(complex<double>(0.0, teig[i]));
     for_each(amat->element_ptr(0,i), amat->element_ptr(0,i+1), [&ex](complex<double>& a) { a *= ex; });
   }
   auto expa = make_shared<ZMatrix>(*amat ^ *amat_sav);

   coeff_ = make_shared<const ZMatrix>(*coeff_ * *expa);

   // print out...
   print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

  }
//
//  // block diagonalize coeff_ in nclosed and nvirt
//  if (nact_)
//    coeff_ = semi_canonical_orb();
//
  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  fci_->update(coeff_);
  fci_->compute();
  fci_->compute_rdm12();

}


void ZSuperCI::one_body_operators(shared_ptr<ZMatrix>& f, shared_ptr<ZMatrix>& fact, shared_ptr<ZMatrix>& factp, shared_ptr<ZMatrix>& gaa,
                                  shared_ptr<ZRotFile>& denom) {
  // calculate 1RDM in an original basis set
  shared_ptr<const ZMatrix> rdm1 = nact_ ? transform_rdm1() : nullptr;
  // make natural orbitals, update coeff_ and transform rdm1
  shared_ptr<ZMatrix> natorb_coeff = make_natural_orbitals(rdm1);
  rdm1 = natorb_rdm1_transform(natorb_coeff, rdm1);

  assert(coeff_->mdim()== nbasis_*2);
  // qvec
  shared_ptr<const ZMatrix> qvec;
  if (nact_) {
    qvec = make_shared<ZQvec>(nbasis_, nact_, geom_, coeff_, nclosed_, fci_, gaunt_, breit_);
    // transform to natural orbitals
    qvec = update_qvec(qvec, natorb_coeff);
  }

  shared_ptr<const ZMatrix> cfock;
  { // Fock operators
    // closed Fock - same as inactive fock
    shared_ptr<const ZMatrix> cfockao = nclosed_ ? make_shared<const DFock>(geom_, hcore_, coeff_->slice(0,nclosed_*2), gaunt_, breit_, /*store half*/false, /*robust*/breit_) : hcore_;
    cfock = make_shared<ZMatrix>(*coeff_ % *cfockao * *coeff_);
    // active Fock operator
    shared_ptr<const ZMatrix> afock;
    if (nact_) {
      shared_ptr<const ZMatrix> afockao = active_fock(rdm1);
      afock = make_shared<ZMatrix>(*coeff_ % *afockao * *coeff_);
    } else {
      afock = make_shared<ZMatrix>(nbasis_*2, nbasis_*2);
    }
    f = make_shared<ZMatrix>(*cfock + *afock);
  }
  { // active-x Fock operator : D_ts cfock_sx + Q_tx
    fact = qvec->copy(); // nbasis_ runs first
    for (int i = 0; i != nact_*2; ++i)
      zaxpy_(nbasis_*2, occup_[i], cfock->element_ptr(0,nclosed_*2+i), 1, fact->data()+i*nbasis_*2, 1);
  }
  { // active Fock' operator (Fts+Fst) / (ns+nt)
    factp = make_shared<ZMatrix>(nact_*2, nact_*2);
    for (int i = 0; i != nact_*2; ++i) {
      for (int j = 0; j != nact_*2; ++j) {
        if (occup_[i] + occup_[j] > zoccup_thresh)
          factp->element(j,i) = (fact->element(j+nclosed_*2,i)+fact->element(i+nclosed_*2,j)) / (occup_[i]+occup_[j]);
        else
          factp->element(j,i) = complex<double> (0.0, 0.0);
      }
    }
  }

  // G matrix (active-active) D_rs,tu Factp_tu - delta_rs nr sum_v Factp_vv
  gaa = factp->clone();
  auto nat_rdm2 = natorb_rdm2_transform(natorb_coeff, fci_->rdm2_av());
  zgemv_("N", nact_*nact_*4, nact_*nact_*4, 1.0, nat_rdm2->data(), nact_*nact_*4, factp->data(), 1, 0.0, gaa->data(), 1);
  complex<double> p = complex<double> (0.0,0.0);
  for (int i = 0; i != nact_*2; ++i) p += occup_[i] * factp->element(i,i);
  for (int i = 0; i != nact_*2; ++i) gaa->element(i,i) -= occup_[i] * p;

  // diagonal denom
  auto dtmp = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);

  complex<double>* target = dtmp->ptr_va();
  for (int i = 0; i != nact_*2; ++i) {
    if (occup_[i] > zoccup_thresh) {
      for (int j = 0; j != nvirt_*2; ++j, ++target)
        *target = (gaa->element(i,i) + occup_[i]*f->element(j+nocc_*2, j+nocc_*2)) / (occup_[i]);
    } else {
      for (int j = 0; j != nvirt_*2; ++j, ++target)
        *target = 1.0/zoccup_thresh;
    } 
  } 

  target = dtmp->ptr_vc();
  for (int i = 0; i != nclosed_*2; ++i)
    for (int j = 0; j != nvirt_*2; ++j, ++target)
      *target = (f->element(j+nocc_*2, j+nocc_*2) - f->element(i, i)) / 1.0; // TODO : check if this factor is 2.0 or 1.0 ?

  target = dtmp->ptr_ca();
  for (int i = 0; i != nact_*2; ++i) {
    if (1.0-occup_[i] > zoccup_thresh) { // TODO : check if the factor is 2.0 - or 1.0 - ...
      for (int j = 0; j != nclosed_*2; ++j, ++target)
        *target = ((f->element(nclosed_*2+i,nclosed_*2+i)-fact->element(i+nclosed_*2,i)) - f->element(j, j)*(1.0-occup_[i])) / (1.0-occup_[i]); // TODO : check on the factors of 2.0
    } else {
      for (int j = 0; j != nclosed_*2; ++j, ++target)
        *target = 1.0/zoccup_thresh;
    }
  }
  const double thresh = 1.0e-8;
  for (int i = 0; i != dtmp->size(); ++i)
    if (fabs(dtmp->data(i)) < thresh) {
      dtmp->data(i) = 1.0e10;
    }

  denom = dtmp;
}

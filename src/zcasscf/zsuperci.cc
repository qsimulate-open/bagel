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


#include <src/rel/dfock.h>
#include <src/zcasscf/zsuperci.h>
#include <src/zcasscf/zsupercimicro.h>
#include <src/rel/reloverlap.h>
#include <src/math/hpw_diis.h>

using namespace std;
using namespace bagel;

void ZSuperCI::compute() {
  shared_ptr<HPW_DIIS<ZMatrix,ZMatrix>> diis;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  double gradient = 1.0e10;

  cout << "     See casscf.log for further information on FCI output " << endl << endl;
  mute_stdcout();
  double orthonorm;
  {
    auto unit = coeff_->clone(); unit->unit();
    orthonorm = ((*coeff_ % *overlap_ * *coeff_) - *unit).rms();
    if (orthonorm > 2.5e-13) throw logic_error("Coefficient is not sufficiently orthnormal.");
  }
  for (int iter = 0; iter != max_iter_; ++iter) {

    if (iter >= diis_start_ && gradient < 1.0e-2 && diis == nullptr) {
      shared_ptr<ZMatrix> tmp = make_shared<ZMatrix>(coeff_->ndim(),coeff_->mdim()/2);
      tmp->copy_block(0, 0, coeff_->ndim(), nocc_*2, coeff_->slice(0, nocc_*2));
      tmp->copy_block(0, nocc_*2, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2,nocc_*2+nvirtnr_));
      tmp->copy_block(0, nocc_*2+nvirtnr_, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2+nvirt_,nocc_*2+nvirt_+nvirtnr_));
      shared_ptr<ZMatrix> unit = make_shared<ZMatrix>(coeff_->mdim()/2, coeff_->mdim()/2);
      unit->unit();
      diis = make_shared<HPW_DIIS<ZMatrix, ZMatrix>>(10, tmp, unit);
    }

    // first perform CASCI to obtain RDMs
    if (nact_) {
      if (iter) fci_->update(coeff_, /*restricted*/true);
      Timer fci_time(0);
      cout << " Executing FCI calculation in Cycle " << iter << endl;
      fci_->compute();
      fci_time.tick_print("ZFCI");
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      fci_time.tick_print("RDMs");
      energy_ = fci_->energy();
    }
    auto grad = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirtnr_*2);

    // compute one-body operators
    shared_ptr<ZMatrix> f, fact, factp, gaa;
    shared_ptr<ZRotFile> denom;
    Timer onebody(0);
    one_body_operators(f, fact, factp, gaa, denom);

    // first, <proj|H|0> is computed
    grad->zero();
    // <a/i|H|0> = f_ai
    grad_vc(f, grad);
    // <a/r|H|0> = cfock_ar n_r + ((as|tu)D_rs,tu)^* = fact_ar
    grad_va(fact, grad);
    // <r/i|H|0> = f_ri - f^inact_is d_sr - (is|tu)P_rs,tu = f_ri - fact_ri
    grad_ca(f, fact, grad);

    onebody.tick_print("One body operators");

    if (!nact_) { // compute energy
      assert(nstate_ == 1 && energy_.size() == 1);
      auto hcoremo = make_shared<ZMatrix>(coeff_->slice(0,nclosed_*2) % *hcore_ * coeff_->slice(0,nclosed_*2));
      *hcoremo += *f->get_submatrix(0, 0, nclosed_*2, nclosed_*2);
      double etmp = 0.0;
      for (int j=0; j!= nclosed_*2; ++j)
        etmp += 0.5 * hcoremo->element(j,j).real();
      etmp += geom_->nuclear_repulsion();
      energy_[0] = etmp;
    }

    // setting error of macro iteration
    gradient = grad->rms();
    if (gradient < thresh_) {
      resume_stdcout();
      // print out...
      if (!nact_) {
        print_iteration(iter, 0, 0, energy_, gradient, timer.tick());
      } else {
        print_iteration(iter, 0, 0, energy_, gradient, timer.tick());
      }
      rms_grad_ = gradient;
      cout << " " << endl;
      cout << "    * Super CI optimization converged. *    " << endl << endl;
      mute_stdcout();
      break;
    }

  // ============================
     // Micro-iterations go here
  // ============================
    shared_ptr<const ZRotFile> cc;
    {
      Timer microiter_time(0);
      ZSuperCIMicro micro(shared_from_this(), grad, denom, f, fact, factp, gaa);
      micro.compute();
      cc = micro.cc();
      microiter_time.tick_print("Microiterations");
    }

    // orbital rotation matrix
    shared_ptr<ZMatrix> amat = cc->unpack<ZMatrix>();
    kramers_adapt(amat, nvirtnr_);
    // multiply -i to make amat hermite (will be compensated), sqrt(2) to recover non-rel limit
    *amat *= sqrt(2.0) * complex<double>(0.0, -1.0);
    VectorB teig(amat->ndim());
    amat->diagonalize(teig);
    auto amat_sav = amat->copy();
    for (int i = 0; i != amat->ndim(); ++i) {
      complex<double> ex = exp(complex<double>(0.0, teig(i)));
      for_each(amat->element_ptr(0,i), amat->element_ptr(0,i+1), [&ex](complex<double>& a) { a *= ex; });
    }
    auto expa = make_shared<ZMatrix>(*amat ^ *amat_sav);
    expa->purify_unitary();

    if (diis == nullptr) {
      auto ctmp = make_shared<ZMatrix>(coeff_->ndim(), nbasis_);
      ctmp->copy_block(0, 0, coeff_->ndim(), nocc_*2, coeff_->slice(0, nocc_*2));
      ctmp->copy_block(0, nocc_*2, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2,nocc_*2+nvirtnr_));
      ctmp->copy_block(0, nocc_*2+nvirtnr_, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2+nvirt_,nocc_*2+nvirt_+nvirtnr_));
      *ctmp *= *expa;
      auto ctmp2 = coeff_->copy();
      ctmp2->copy_block(0, 0, coeff_->ndim(), nocc_*2, ctmp->slice(0, nocc_*2));
      ctmp2->copy_block(0, nocc_*2, coeff_->ndim(), nvirtnr_, ctmp->slice(nocc_*2,nocc_*2+nvirtnr_));
      ctmp2->copy_block(0, nocc_*2+nvirt_, coeff_->ndim(), nvirtnr_, ctmp->slice(nocc_*2+nvirtnr_,nocc_*2+nvirtnr_*2));
      coeff_ = make_shared<const ZMatrix>(*ctmp2);
    } else {
      shared_ptr<const ZMatrix> mcc = diis->extrapolate(expa);
      auto ctmp2 = coeff_->copy();
      ctmp2->copy_block(0, 0, coeff_->ndim(), nocc_*2, mcc->slice(0, nocc_*2));
      ctmp2->copy_block(0, nocc_*2, coeff_->ndim(), nvirtnr_, mcc->slice(nocc_*2,nocc_*2+nvirtnr_));
      ctmp2->copy_block(0, nocc_*2+nvirt_, coeff_->ndim(), nvirtnr_, mcc->slice(nocc_*2+nvirtnr_,nocc_*2+nvirtnr_*2));
      coeff_ = make_shared<const ZMatrix>(*ctmp2);
    }

    // synchronization
    mpi__->broadcast(const_pointer_cast<ZMatrix>(coeff_)->data(), coeff_->size(), 0);

    // print out...
    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());
    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      if (real(rms_grad_) > thresh_) cout << "    * The calculation did NOT converge. *    " << endl;
      cout << "    * Max iteration reached in the Super CI macro interations. *     " << endl << endl;
    }
    {
      auto unit = coeff_->clone(); unit->unit();
      auto orthonorm2 = ((*coeff_ % *overlap_ * *coeff_) - *unit).rms();
      if (orthonorm2 / orthonorm > 1.0e+01)
        throw logic_error("should not happen");
    }
    mute_stdcout();

  }
  resume_stdcout();

  // TODO : block diagonalize coeff_ in nclosed and nvirt

  // the following is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  if (nact_) {
    fci_->update(coeff_, /*restricted*/true);
    fci_->compute();
    fci_->compute_rdm12();
  }

}

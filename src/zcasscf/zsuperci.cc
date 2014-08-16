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
      Timer fci_time(0);
      if (iter) fci_->update(coeff_, /*restricted*/true);
      cout << " Executing FCI calculation in Cycle " << iter << endl;
      fci_->compute();
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      energy_ = fci_->energy();
      fci_time.tick_print("FCI and RDMs");
    }
    auto grad = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirtnr_*2);

    // compute one-boedy operators
    shared_ptr<ZMatrix> f, fact, factp, gaa;
    shared_ptr<ZRotFile> denom;
    Timer onebody(0);
    one_body_operators(f, fact, factp, gaa, denom);
    onebody.tick_print("One body operators");

    // first, <proj|H|0> is computed
    grad->zero();
    // <a/i|H|0> = f_ai
    grad_vc(f, grad);
    // <a/r|H|0> = cfock_ar n_r + ((as|tu)D_rs,tu)^* = fact_ar
    grad_va(fact, grad);
    // <r/i|H|0> = f_ri - f^inact_is d_sr - (is|tu)P_rs,tu = f_ri - fact_ri
    grad_ca(f, fact, grad);

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
      cout << "    * ZSuperCI optimization converged    " << endl << endl;
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

    // print out...
    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());
    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      if (real(rms_grad_) > thresh_) cout << "    * The calculation did NOT converge. *    " << endl;
      cout << "    * Max iteration reached in the Super CI macro interations. *     " << endl << endl;
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


void ZSuperCI::one_body_operators(shared_ptr<ZMatrix>& f, shared_ptr<ZMatrix>& fact, shared_ptr<ZMatrix>& factp, shared_ptr<ZMatrix>& gaa,
                                  shared_ptr<ZRotFile>& denom) {
  bool a2approx = idata_->get<bool>("a2approx", false);
  assert(coeff_->mdim()== nbasis_*2);

  // qvec ; electronic contributions only
  shared_ptr<const ZMatrix> qvec;
  if (nact_) {
    // extract electronic orbitals from coeff
    auto coefftmp = make_shared<ZMatrix>(coeff_->ndim(), nbasis_);
    coefftmp->copy_block(0, 0, coeff_->ndim(), nocc_*2, coeff_->slice(0, nocc_*2));
    coefftmp->copy_block(0, nocc_*2, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2, nocc_*2+nvirtnr_));
    coefftmp->copy_block(0, nocc_*2+nvirtnr_, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2+nvirt_, nocc_*2+nvirt_+nvirtnr_));
    if (!a2approx) {
      qvec = make_shared<ZQvec>(nbasis_, nact_, geom_, coefftmp, coefftmp->slice_copy(nclosed_*2, nocc_*2), nclosed_, fci_, gaunt_, breit_);
    } else {
      qvec = make_shared<const ZMatrix>(nocc_*2+nvirtnr_*2, nact_*2);
    }
  }

  // calculate 1RDM in an original basis set
  shared_ptr<const ZMatrix> rdm1 = nact_ ? transform_rdm1() : nullptr;
  // make natural orbitals, update coeff_ and transform rdm1
  shared_ptr<ZMatrix> natorb_coeff;
  if (nact_) {
    natorb_coeff = make_natural_orbitals(rdm1);
    coeff_ = update_coeff(coeff_, natorb_coeff);
    qvec = update_qvec(qvec, natorb_coeff);
    rdm1 = natorb_rdm1_transform(natorb_coeff, rdm1);
  }

  shared_ptr<const ZMatrix> cfock;
  { // Fock operators
    // extract electronic orbitals from coeff
    auto coefftmp = make_shared<ZMatrix>(coeff_->ndim(), nbasis_);
    coefftmp->copy_block(0, 0, coeff_->ndim(), nocc_*2, coeff_->slice(0, nocc_*2));
    coefftmp->copy_block(0, nocc_*2, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2, nocc_*2+nvirtnr_));
    coefftmp->copy_block(0, nocc_*2+nvirtnr_, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2+nvirt_, nocc_*2+nvirt_+nvirtnr_));

    // closed Fock - same as inactive fock
    if (!nact_) {
      shared_ptr<const ZMatrix> cfockao = nclosed_ ? make_shared<const DFock>(geom_, hcore_, coeff_->slice_copy(0,nclosed_*2), gaunt_, breit_, /*store half*/false, /*robust*/breit_) : hcore_;
      cfock = make_shared<ZMatrix>(*coefftmp % *cfockao * *coefftmp);
    } else {
      cfock = make_shared<const ZMatrix>(*coefftmp % *fci_->jop()->core_fock() * *coefftmp);
    }
    // active Fock operator
    shared_ptr<const ZMatrix> afock;
    if (nact_) {
      shared_ptr<const ZMatrix> afockao = active_fock(rdm1);
      afock = make_shared<ZMatrix>(*coefftmp % *afockao * *coefftmp);
    } else {
      afock = cfock->clone();
    }
    f = make_shared<ZMatrix>(*cfock + *afock);
  }
  if (nact_) { // x-active Fock operator : cfock_xs^ n_s + Q_xt^*
    fact = qvec->get_conjg();
    for (int i = 0; i != nact_*2; ++i)
      zaxpy_(qvec->ndim(), occup_[i], cfock->element_ptr(0,nclosed_*2+i), 1, fact->data()+i*qvec->ndim(), 1);
  }
  if (nact_) { // active Fock' operator (Fts+Fst) / (ns+nt)
    factp = make_shared<ZMatrix>(nact_*2, nact_*2);
    shared_ptr<ZMatrix> fact_conjg = fact->get_conjg();
    for (int i = 0; i != nact_*2; ++i) {
      for (int j = 0; j != nact_*2; ++j) {
        if (occup_[i] + occup_[j] > zoccup_thresh)
          factp->element(j,i) = (fact->element(j+nclosed_*2,i)+fact_conjg->element(i+nclosed_*2,j)) / (occup_[i]+occup_[j]);
        else
          factp->element(j,i) = complex<double> (0.0, 0.0);
      }
    }
  }

  // G matrix (active-active) D_rs,tu Factp_tu - delta_rs nr sum_v Factp_vv
  if (nact_) {
    gaa = factp->clone();
    shared_ptr<const ZMatrix> nat_rdm2 = natorb_rdm2_transform(natorb_coeff, fci_->rdm2_av());
    if (!a2approx)
      zgemv_("N", nact_*nact_*4, nact_*nact_*4, 1.0, nat_rdm2->data(), nact_*nact_*4, factp->get_conjg()->data(), 1, 0.0, gaa->data(), 1);
    complex<double> p = complex<double> (0.0,0.0);
    for (int i = 0; i != nact_*2; ++i) p += occup_[i] * factp->element(i,i);
    for (int i = 0; i != nact_*2; ++i) gaa->element(i,i) -= occup_[i] * p;
  }

  // diagonal denom
  {
    int nvirt_tmp = nvirtnr_;
    auto dtmp = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_tmp*2);

    complex<double>* target = dtmp->ptr_va();
    for (int i = 0; i != nact_*2; ++i) {
      if (occup_[i] > zoccup_thresh) {
        for (int j = 0; j != nvirt_tmp*2; ++j, ++target)
          *target = (gaa->element(i,i) + occup_[i]*f->element(j+nocc_*2, j+nocc_*2)) / (occup_[i]);
      } else {
        for (int j = 0; j != nvirt_tmp*2; ++j, ++target)
          *target = 1.0/zoccup_thresh;
      }
    }

    target = dtmp->ptr_vc();
    for (int i = 0; i != nclosed_*2; ++i)
      for (int j = 0; j != nvirt_tmp*2; ++j, ++target)
        *target = (f->element(j+nocc_*2, j+nocc_*2) - f->element(i, i)) / 2.0; // 2.0 to recover non-rel limit

    if (nact_) {
      target = dtmp->ptr_ca();
      for (int i = 0; i != nact_*2; ++i) {
        if (1.0-occup_[i] > zoccup_thresh) {
          for (int j = 0; j != nclosed_*2; ++j, ++target)
            *target = ((f->element(nclosed_*2+i,nclosed_*2+i)-fact->element(i+nclosed_*2,i)) - f->element(j, j)*(1.0-occup_[i])) / (1.0-occup_[i]);
        } else {
          for (int j = 0; j != nclosed_*2; ++j, ++target)
            *target = 1.0/zoccup_thresh;
        }
      }
    }
    const double thresh = 1.0e-8;
    for (int i = 0; i != dtmp->size(); ++i)
      if (fabs(dtmp->data(i)) < thresh) {
        dtmp->data(i) = 1.0e10;
      }

    denom = dtmp;
  }
}

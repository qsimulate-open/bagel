//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zuperci.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson Bates <jefferson.bates@northwestern.edu>
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
#include <src/mat1e/rel/reloverlap.h>
#include <src/multi/zcasscf/zqvec.h>
#include <src/multi/zcasscf/zsuperci.h>
#include <src/multi/zcasscf/zsupercimicro.h>
#include <src/util/math/hpw_diis.h>
#include <src/util/math/step_restrict_bfgs.h>
#include <src/scf/dhf/population_analysis.h>

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
  muffle_->mute();
  double orthonorm;
  {
    auto unit = coeff_->clone(); unit->unit();
    orthonorm = ((*coeff_ % *overlap_ * *coeff_) - *unit).rms();
    if (orthonorm > 2.5e-13)
      cout << "Coefficient is not sufficiently orthnormal: " << setprecision(10) << setw(15) << orthonorm << endl;;
  }
  for (int iter = 0; iter != max_iter_; ++iter) {

    // DIIS setup
    if (iter >= diis_start_ && gradient < 1.0e-2 && diis == nullptr) {
      shared_ptr<ZMatrix> unit = make_shared<ZMatrix>(coeff_->mdim()/2, coeff_->mdim()/2);
      unit->unit();
      diis = make_shared<HPW_DIIS<ZMatrix, ZMatrix>>(10, coeff_->electronic_part(), unit);
    }

    // first perform CASCI to obtain RDMs
    if (nact_) {
      if (iter) fci_->update(coeff_);
      Timer fci_time(0);
      cout << " Executing FCI calculation in Cycle " << iter << endl;
      fci_->compute();
      fci_time.tick_print("ZFCI");
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      fci_time.tick_print("RDMs");
      if (iter > 0) prev_energy_ = energy_;
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
      muffle_->unmute();
      // print out...
      print_iteration(iter, energy_, gradient, timer.tick());
      rms_grad_ = gradient;
      cout << endl;
      // output energy change for last cycle
      cout << "    * State averaged energy change from last cycle = " << setprecision(6) << scientific
           << blas::average(prev_energy_) - blas::average(energy_) << endl;
      cout << "    * Super CI optimization converged. *    " << endl << endl;
      muffle_->mute();
      break;
    }

  // ============================
  //   Micro-iterations go here
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
    shared_ptr<ZMatrix> amat = cc->unpack();
    if (tsymm_)
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
      // rotate electronic orbitals
      shared_ptr<ZMatrix> ctmp = coeff_->electronic_part();
      *ctmp = *ctmp * *expa;
      coeff_ = coeff_->update_electronic(ctmp);
    } else {
      // DIIS extrapolate
      shared_ptr<const ZMatrix> mcc = diis->extrapolate(expa);
      coeff_ = coeff_->update_electronic(mcc);
    }

    // synchronization
    mpi__->broadcast(const_pointer_cast<RelCoeff_Block>(coeff_)->data(), coeff_->size(), 0);

    // print out...
    muffle_->unmute();
    print_iteration(iter, energy_, gradient, timer.tick());
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
    muffle_->mute();

  }
  muffle_->unmute();

  // TODO : block diagonalize coeff_ in nclosed and nvirt

  // the following is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  if (nact_) {
    fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
  }

  // print out orbital populations, if needed
  if (idata_->get<bool>("pop", false)) {
    Timer pop_timer;
    cout << " " << endl;
    cout << "    * Printing out population analysis of super-CI optimized orbitals to casscf.log" << endl;
    muffle_->mute();
    population_analysis(geom_, coeff_->striped_format()->slice(0, 2*(nclosed_+nact_+nvirtnr_)), overlap_, tsymm_);
    muffle_->unmute();
    pop_timer.tick_print("population analysis");
  }
}


void ZSuperCI::one_body_operators(shared_ptr<ZMatrix>& f, shared_ptr<ZMatrix>& fact, shared_ptr<ZMatrix>& factp, shared_ptr<ZMatrix>& gaa,
                                  shared_ptr<ZRotFile>& denom) {
  assert(coeff_->mdim()== nbasis_*2);

  // qvec ; electronic contributions only
  shared_ptr<const ZMatrix> qvec;
  if (nact_) {
    qvec = make_shared<ZQvec>(nbasis_, nact_, geom_, coeff_->electronic_part(), coeff_->slice_copy(nclosed_*2, nocc_*2), nclosed_, fci_, gaunt_, breit_);
  }

  // calculate 1RDM in an original basis set
  shared_ptr<const ZMatrix> rdm1 = nact_ ? fci_->rdm1_av() : nullptr;
  // make natural orbitals, update coeff_ and transform rdm1
  shared_ptr<ZMatrix> natorb_coeff;
  if (nact_) {
    pair<shared_ptr<ZMatrix>, VectorB> natorb_tmp = make_natural_orbitals(rdm1);
    occup_ = natorb_tmp.second;
    natorb_coeff = natorb_tmp.first;
    coeff_ = update_coeff(coeff_, natorb_coeff);
    qvec = update_qvec(qvec, natorb_coeff);
    rdm1 = make_shared<ZMatrix>(*natorb_coeff % *rdm1 * *natorb_coeff);
  }

  shared_ptr<const ZMatrix> cfock;
  { // Fock operators
    shared_ptr<const ZMatrix> coeff_elec = coeff_->electronic_part();

    // closed Fock - same as inactive fock
    if (!nact_) {
      shared_ptr<const ZMatrix> cfockao = nclosed_ ? make_shared<const DFock>(geom_, hcore_, coeff_->slice_copy(0,nclosed_*2), gaunt_, breit_, /*store half*/false, /*robust*/breit_) : hcore_;
      cfock = make_shared<ZMatrix>(*coeff_elec % *cfockao * *coeff_elec);
    } else {
      cfock = make_shared<const ZMatrix>(*coeff_elec % *fci_->jop()->core_fock() * *coeff_elec);
    }
    // active Fock operator
    shared_ptr<const ZMatrix> afock;
    if (nact_) {
      shared_ptr<const ZMatrix> afockao = compute_active_fock(coeff_->slice(nclosed_*2, nocc_*2), rdm1);
      afock = make_shared<ZMatrix>(*coeff_elec % *afockao * *coeff_elec);
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

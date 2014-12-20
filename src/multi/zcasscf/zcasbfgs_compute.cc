//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf_compute.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/multi/zcasscf/zqvec.h>
#include <src/multi/zcasscf/zcasbfgs.h>
#include <src/util/math/step_restrict_bfgs.h>
#include <src/scf/dhf/dfock.h>
#include <src/mat1e/rel/reloverlap.h>

using namespace std;
using namespace bagel;


void ZCASBFGS::compute() {
  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.
  shared_ptr<SRBFGS<ZRotFile>> ele_srbfgs;
  shared_ptr<SRBFGS<ZRotFile>> pos_srbfgs;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  // allocate unitary rotation matrices
  auto ele_x = make_shared<ZMatrix>((nocc_+nvirtnr_)*2, (nocc_+nvirtnr_)*2);
  ele_x->unit();
  vector<double> ele_energy;
  auto pos_x = make_shared<ZMatrix>(nocc_*2 + nneg_, nocc_*2 + nneg_);
  pos_x->unit();
  vector<double> pos_energy;
  bool ele_conv = false;
  bool pos_conv = false;

//  auto cold = coeff_->clone(); // TODO : needed if step rejection is implemented in the future

  bool optimize_electrons = idata_->get<bool>("optimize_electrons", true);
  const bool only_electrons = idata_->get<bool>("only_electrons", false);
  if (only_electrons)  cout << "     Orbital optimization for electronic orbitals only " << endl;
  double orthonorm;
  {
    auto unit = coeff_->clone(); unit->unit();
    orthonorm = ((*coeff_ % *overlap_ * *coeff_) - *unit).rms();
    if (orthonorm > 2.5e-13) throw logic_error("Coefficient is not sufficiently orthnormal.");
  }
  cout << "     See casscf.log for further information on FCI output " << endl << endl;
  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    if (nact_) {
      if (iter) fci_->update(coeff_, /*restricted*/true);
      cout << " Executing FCI calculation in Cycle " << iter << endl;
      Timer fci_time(0);
      fci_->compute();
      fci_time.tick_print("ZFCI");
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      fci_time.tick_print("RDMs");
    }

    // TODO : compute one body operators only for the subspace being optimized; presently full coefficient is used to transform from AO to MO
    // TODO : make one body operators function to simplify the driver ; could combine with super-CI, but will take some work
    Timer onebody(0);
    // calculate 1RDM in an original basis set
    shared_ptr<const ZMatrix> rdm1 = nact_ ? transform_rdm1() : nullptr;

    // closed Fock operator
    shared_ptr<const ZMatrix> cfock;
    shared_ptr<const ZMatrix> cfockao;
    if (!nact_) {
      cfockao = nclosed_ ? make_shared<const DFock>(geom_, hcore_, coeff_->slice_copy(0,nclosed_*2), gaunt_, breit_, /*store half*/false, /*robust*/breit_) : hcore_;
      cfock = make_shared<ZMatrix>(*coeff_ % *cfockao * *coeff_);
    } else {
      cfockao = fci_->jop()->core_fock();
      cfock = make_shared<const ZMatrix>(*coeff_ % *cfockao * *coeff_);
    }

    // active Fock operator
    shared_ptr<const ZMatrix> afock;
    if (nact_) {
      shared_ptr<const ZMatrix> afockao = active_fock(rdm1, /*with_hcore*/false, /*bfgs*/true);
      afock = make_shared<ZMatrix>(*coeff_ % *afockao * *coeff_);
    } else {
      afock = make_shared<ZMatrix>(nbasis_*2, nbasis_*2);
    }
    assert(coeff_->mdim()== nbasis_*2);

    // qvec
    shared_ptr<const ZMatrix> qvec;
    if (nact_) {
      qvec = make_shared<ZQvec>(nbasis_*2, nact_, geom_, coeff_, coeff_->slice_copy(nclosed_*2,nocc_*2), nclosed_, fci_, gaunt_, breit_)->get_conjg();
    }

    // compute orbital gradients
    shared_ptr<ZRotFile> grad = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);
    grad_vc(cfock, afock, grad);
    grad_va(cfock, qvec, rdm1, grad);
    grad_ca(cfock, afock, qvec, rdm1, grad);
    *grad *= 2.0;
    if (pos_conv) {
      auto newgrad = copy_positronic_rotations(grad);
      cout << setprecision(4) << scientific << " Positron gradient RMS = " << newgrad->rms() << endl;
      if (newgrad->rms() > thresh_) pos_conv = false;
      if (newgrad->rms() > thresh_) cout << " POSITRONS NOT CONVERGED " << endl;
    }

    // compute approximate diagonal hessian
    if (iter == 0) {
      shared_ptr<ZRotFile> denom = compute_denom(cfock, afock, qvec, rdm1);
      // IMPROVISED LEVEL SHIFT
      const bool shift = idata_->get<bool>("shift", false);
      if (shift) {
        level_shift_ = find_level_shift(denom);
        shared_ptr<ZRotFile> diagonal_shift = denom->clone();
        diagonal_shift->fill(level_shift_);
        denom->ax_plus_y(-1.0, diagonal_shift);
      }
      { // electronic rotation bfgs
        auto newdenom = copy_electronic_rotations(denom);
        ele_srbfgs = make_shared<SRBFGS<ZRotFile>>(newdenom);
      }
      { // positronic rotation bfgs
        auto newdenom = copy_positronic_rotations(denom);
        const double thresh = 1.0e-8;
        for (int i = 0; i != newdenom->size(); ++i)
          if (fabs(newdenom->data(i)) < thresh) {
            newdenom->data(i) = 1.0e10;
          }
        pos_srbfgs = make_shared<SRBFGS<ZRotFile>>(newdenom);
      }
    }
    onebody.tick_print("One body operators");

    // compute unitary rotation matrix for given subspace, extract subspace gradient, and update energy vectors
    shared_ptr<ZMatrix> expa;
    if (optimize_electrons) {
      expa = compute_unitary_rotation(ele_energy, ele_srbfgs, ele_x, nvirtnr_, cfockao, grad, optimize_electrons);
    } else {
      expa = compute_unitary_rotation(pos_energy, pos_srbfgs, pos_x, nneg_/2,  cfockao, grad, optimize_electrons);
    }

    const double gradient = grad->rms();

//    cold = coeff_->copy(); // TODO : copy old coefficient if step rejection is ever implemented
    if (optimize_electrons) {
      // extract electronic orbitals from coefficient
      auto ctmp = make_shared<ZMatrix>(coeff_->ndim(), coeff_->mdim()/2);
      ctmp->copy_block(0, 0, coeff_->ndim(), nocc_*2 + nvirtnr_, coeff_->slice(0, nocc_*2 + nvirtnr_));
      ctmp->copy_block(0, nocc_*2 + nvirtnr_, coeff_->ndim(), nvirtnr_, coeff_->slice(nocc_*2 + nvirt_, nocc_*2 + nvirt_ + nvirtnr_));
      // rotate orbitals
      *ctmp = *ctmp * *expa;
      // copy back to full coeff
      auto ctmp2 = coeff_->copy();
      ctmp2->copy_block(0, 0, coeff_->ndim(), nocc_*2 + nvirtnr_, ctmp->slice(0, nocc_*2 + nvirtnr_));
      ctmp2->copy_block(0, nocc_*2 + nvirt_, coeff_->ndim(), nvirtnr_, ctmp->slice(nocc_*2 + nvirtnr_, ctmp->mdim()));
      coeff_ = make_shared<const ZMatrix>(*ctmp2);
    } else {
      // extract occupied and positronic orbitals from coefficient
      auto ctmp = make_shared<ZMatrix>(coeff_->ndim(), coeff_->mdim()/2 + nocc_*2);
      ctmp->copy_block(0, 0, coeff_->ndim(), nocc_*2, coeff_->slice(0, nocc_*2));
      ctmp->copy_block(0, nocc_*2, coeff_->ndim(), nneg_/2, coeff_->slice(nocc_*2 + nvirtnr_, nocc_*2 + nvirt_));
      ctmp->copy_block(0, nocc_*2 + nneg_/2, coeff_->ndim(), nneg_/2, coeff_->slice(nocc_*2 + nvirt_ + nvirtnr_, nocc_*2 + nvirt_*2));
      // rotate orbitals
      *ctmp = *ctmp * *expa;
      // copy back to full coeff
      auto ctmp2 = coeff_->copy();
      ctmp2->copy_block(0, 0, coeff_->ndim(), nocc_*2, ctmp->slice(0, nocc_*2));
      ctmp2->copy_block(0, nocc_*2 + nvirtnr_, coeff_->ndim(), nneg_/2, ctmp->slice(nocc_*2, nocc_*2 +nneg_/2));
      ctmp2->copy_block(0, nocc_*2 + nvirtnr_ + nvirt_, coeff_->ndim(), nneg_/2, ctmp->slice(nocc_*2 + nneg_/2, ctmp->mdim()));
      coeff_ = make_shared<const ZMatrix>(*ctmp2);
    }
    // for next BFGS extrapolation
    if (optimize_electrons) {
      *ele_x *= *expa;
    } else {
      *pos_x *= *expa;
    }

    // synchronization
    mpi__->broadcast(const_pointer_cast<ZMatrix>(coeff_)->data(), coeff_->size(), 0);

    // print energy
    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    // Set logic flags based upon convergence criteria and switch optimization subspaces accordingly
    if (!optimize_electrons) {
      // end of e-p iteration
      if (gradient < thresh_) pos_conv = true; // positrons converged
      optimize_electrons = ele_conv ? false : true; // switch to electrons if e-e rotations are NOT converged
    } else {
      // end of e-e iteration
      if (gradient < thresh_) ele_conv = true; // electrons converged
      // don't switch to positrons if doing only e-e rotations or positrons are converged
      optimize_electrons = (only_electrons || pos_conv) ? true : false;
    }
    // check convergence
    if ((ele_conv && only_electrons) || (pos_conv && ele_conv)) {
      cout << " " << endl;
      // output energy change for last cycle
      cout << "    * State averaged energy change from last cycle = " << setprecision(6) << scientific
           << blas::average(prev_energy_) - blas::average(energy_) << endl;
      cout << "    * quasi-Newton optimization converged. *   " << endl << endl;
      rms_grad_ = gradient;
      mute_stdcout();
      break;
    }
    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      if (real(rms_grad_) > thresh_) cout << "    * The calculation did NOT converge. *    " << endl;
      cout << "    * Max iteration reached during the quasi-Newton optimization. *     " << endl << endl;
    }
    {
      auto unit = coeff_->clone(); unit->unit();
      auto orthonorm2 = ((*coeff_ % *overlap_ * *coeff_) - *unit).rms();
      if (orthonorm2 / orthonorm > 1.0e+01)
        throw logic_error("Coefficient has lost orthonormality during optimization");
    }
    mute_stdcout();
  }
  if (energy_.size() == 0)
    optimize_electrons == true ? energy_.push_back(ele_energy.back()) : energy_.push_back(pos_energy.back());

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  resume_stdcout();
  if (nact_) {
    fci_->update(coeff_, /*restricted*/true);
    fci_->compute();
    fci_->compute_rdm12();
  }

}

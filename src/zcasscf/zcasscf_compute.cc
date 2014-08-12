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

#include <src/zcasscf/zqvec.h>
#include <src/math/step_restrict_bfgs.h>
#include <src/rel/dfock.h>
#include <src/zcasscf/zcasbfgs.h>
#include <src/rel/reloverlap.h>

#define BOTHSPACES

using namespace std;
using namespace bagel;


void ZCASBFGS::compute() {
  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.
  shared_ptr<SRBFGS<ZRotFile>> srbfgs;
  shared_ptr<SRBFGS<ZRotFile>> ele_srbfgs;
  shared_ptr<SRBFGS<ZRotFile>> pos_srbfgs;
#ifdef BOTHSPACES
  const bool tight = idata_->get<bool>("tight", false);
  const int limmem = idata_->get<int>("limited_memory", 0);
#endif

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

#ifdef BOTHSPACES
  auto x = make_shared<ZMatrix>(nbasis_*2, nbasis_*2);
  x->unit();
#else
  auto ele_x = make_shared<ZMatrix>((nocc_+nvirtnr_)*2, (nocc_+nvirtnr_)*2);
  ele_x->unit();
  vector<double> ele_energy;
  auto pos_x = make_shared<ZMatrix>(nocc_*2 + nneg_, nocc_*2 + nneg_);
  pos_x->unit();
  vector<double> pos_energy;
  bool ele_conv = false;
  bool pos_conv = false;
#endif

  // TODO for debug, we may rotate coefficients. The magnitude can be specified in the input
  const bool ___debug___break_kramers = false;
  const bool ___debug___with_kramers = idata_->get<bool>("debugkramers", true);
  if (___debug___break_kramers)
    ___debug___orbital_rotation(___debug___with_kramers);

  auto cold = coeff_->clone();

  bool optimize_electrons = idata_->get<bool>("optimize_electrons", true);
  cout << " optimizing the electrons ? " << optimize_electrons << endl;
  cout << " See casscf.log for further information on FCI output " << endl;
  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    if (nact_) {
      mute_stdcout(/*fci*/true);
      if (iter) fci_->update(coeff_, /*restricted*/true);
      cout << " Executing FCI calculation in Cycle " << iter << endl;
      fci_->compute();
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      resume_stdcout();
    }

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
      shared_ptr<const ZMatrix> afockao = active_fock(rdm1);
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

    // get energy
    if (nact_) {
#ifdef BOTHSPACES
      energy_.push_back((fci_->energy())[0]);
#else
      optimize_electrons == true ? ele_energy.push_back((fci_->energy())[0]) : pos_energy.push_back((fci_->energy())[0]);
#endif
    } else {
      assert(nstate_ == 1);
#ifdef BOTHSPACES
      energy_.resize(iter+1);
      energy_[iter] = geom_->nuclear_repulsion();
      auto mo = make_shared<ZMatrix>(*coeff_ % (*cfockao+*hcore_) * *coeff_);
      for (int i = 0; i != nclosed_*2; ++i)
        energy_[iter] += 0.5*mo->element(i,i).real();
#else
      optimize_electrons == true ? ele_energy.resize(iter/2+1) : pos_energy.resize((iter-1)/2+1);
      optimize_electrons == true ? ele_energy[iter/2] = geom_->nuclear_repulsion() : pos_energy[(iter-1)/2] = geom_->nuclear_repulsion();
      auto mo = make_shared<ZMatrix>(*coeff_ % (*cfockao+*hcore_) * *coeff_);
      for (int i = 0; i != nclosed_*2; ++i)
        optimize_electrons == true ? ele_energy[iter/2] += 0.5*mo->element(i,i).real() : pos_energy[(iter-1)/2] += 0.5*mo->element(i,i).real();
#endif
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
      srbfgs = make_shared<SRBFGS<ZRotFile>>(denom);
      { // electronic rotation bfgs
        auto newdenom = ___debug___copy_electronic_rotations(denom);
        ele_srbfgs = make_shared<SRBFGS<ZRotFile>>(newdenom);
      }
      { // positronic rotation bfgs
        auto newdenom = ___debug___copy_positronic_rotations(denom);
        const double thresh = 1.0e-8;
        for (int i = 0; i != newdenom->size(); ++i)
          if (fabs(newdenom->data(i)) < thresh) {
            newdenom->data(i) = 1.0e10;
          }
        pos_srbfgs = make_shared<SRBFGS<ZRotFile>>(newdenom);
      }
    }

    // compute orbital gradients
    shared_ptr<ZRotFile> grad = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);
    grad_vc(cfock, afock, grad);
    grad_va(cfock, qvec, rdm1, grad);
    grad_ca(cfock, afock, qvec, rdm1, grad);
    *grad *= 2.0;
    if (!___debug___break_kramers)
      kramers_adapt(grad, nclosed_, nact_, nvirt_);

    if (___debug___break_kramers) {
      ___debug___print_gradient(grad, ___debug___with_kramers);
      ___debug___compute_hessian(cfock, afock, qvec, ___debug___with_kramers);
    }

#ifdef BOTHSPACES
    auto xlog = make_shared<ZRotFile>(x->log(4), nclosed_*2, nact_*2, nvirt_*2);
#else
    shared_ptr<ZRotFile> xlog;
    shared_ptr<ZRotFile> ele_rot;
    shared_ptr<ZRotFile> pos_rot;
    bool reset;
    if (optimize_electrons) {
      cout << "+++ Optimizing electrons +++ " << endl;
      xlog    = make_shared<ZRotFile>(ele_x->log(4), nclosed_*2, nact_*2, nvirtnr_*2);
      tie(ele_rot, ele_energy, grad, xlog, reset) = ___debug___optimize_subspace_rotations(ele_energy, grad, xlog, ele_srbfgs, cold, optimize_electrons);
      kramers_adapt(ele_rot, nclosed_, nact_, nvirtnr_);
    } else {
      cout << "+++ Optimizing positrons +++ " << endl;
      xlog    = make_shared<ZRotFile>(pos_x->log(4), nclosed_*2, nact_*2, nneg_);
      tie(pos_rot, pos_energy, grad, xlog, reset) = ___debug___optimize_subspace_rotations(pos_energy, grad, xlog, pos_srbfgs, cold, optimize_electrons);
      kramers_adapt(pos_rot, nclosed_, nact_, nneg_/2);
    }
#endif

#ifdef BOTHSPACES
    mute_stdcout(/*fci*/false);
    cout << " " << endl;
    cout << " ++++++++++++++++++++++++ " << endl;
    cout << " Starting microiterations " << endl;
    cout << " ++++++++++++++++++++++++ " << endl;
    cout << setprecision(6) << " gradient norm      = " << grad->norm() << endl;
    cout << " " << endl;
    auto reset = srbfgs->check_step(energy_, grad, xlog, tight, limmem);
    if (reset) {
      cout << " STEP DOES NOT MEET PROPER CRITERIA : Please backtrack. " << endl;
    }
    shared_ptr<ZRotFile> a = srbfgs->more_sorensen_extrapolate(grad, xlog);
    resume_stdcout();
    if (!___debug___break_kramers)
      kramers_adapt(a, nclosed_, nact_, nvirt_);
#endif
#ifdef BOTHSPACES
    shared_ptr<ZMatrix> amat = a->unpack<ZMatrix>();
#else
    shared_ptr<ZMatrix> amat;
    if (optimize_electrons) {
      amat = ele_rot->unpack<ZMatrix>();
    } else {
      amat = pos_rot->unpack<ZMatrix>();
    }
#endif

#ifdef BOTHSPACES
    const double gradient = grad->rms();
#else
    const double gradient = grad->rms();//optimize_electrons ? ___debug___copy_electronic_rotations(grad)->rms() : ___debug___copy_positronic_rotations(grad)->rms();
#endif

    // multiply -1 from the formula taken care of in extrap. multiply -i to make amat hermite (will be compensated)
    *amat *= 1.0 * complex<double>(0.0, -1.0);

    // restore the matrix from RotFile
    VectorB teig(amat->ndim());
    amat->diagonalize(teig);
    auto amat_sav = amat->copy();
    for (int i = 0; i != amat->ndim(); ++i) {
      complex<double> ex = exp(complex<double>(0.0, teig(i)));
      for_each(amat->element_ptr(0,i), amat->element_ptr(0,i+1), [&ex](complex<double>& a) { a *= ex; });
    }
    auto expa = make_shared<ZMatrix>(*amat ^ *amat_sav);
    if (!___debug___break_kramers) {
#ifdef BOTHSPACES
      kramers_adapt(expa, nvirt_);
#else
      if (optimize_electrons) {
        kramers_adapt(expa, nvirtnr_);
      } else {
        kramers_adapt(expa, nneg_/2);
      }
#endif
    }

    cold = coeff_->copy();
#ifdef BOTHSPACES
    coeff_ = make_shared<const ZMatrix>(*coeff_ * *expa);
#else
    if (optimize_electrons) {
      int nvirtnr = nvirt_ - nneg_/2;
      auto ctmp = make_shared<ZMatrix>(coeff_->ndim(), coeff_->mdim()/2);
      ctmp->copy_block(0, 0, coeff_->ndim(), nocc_*2 + nvirtnr, coeff_->slice(0, nocc_*2 + nvirtnr));
      ctmp->copy_block(0, nocc_*2 + nvirtnr, coeff_->ndim(), nvirtnr, coeff_->slice(nocc_*2 + nvirt_, nocc_*2 + nvirt_ + nvirtnr));
      *ctmp = *ctmp * *expa;
      auto ctmp2 = coeff_->copy();
      ctmp2->copy_block(0, 0, coeff_->ndim(), nocc_*2 + nvirtnr, ctmp->slice(0, nocc_*2 + nvirtnr));
      ctmp2->copy_block(0, nocc_*2 + nvirt_, coeff_->ndim(), nvirtnr, ctmp->slice(nocc_*2 + nvirtnr, ctmp->mdim()));
      coeff_ = make_shared<const ZMatrix>(*ctmp2);
    } else {
      int nvirtnr = nvirt_ - nneg_/2;
      auto ctmp = make_shared<ZMatrix>(coeff_->ndim(), coeff_->mdim()/2 + nocc_*2);
      ctmp->copy_block(0, 0, coeff_->ndim(), nocc_*2, coeff_->slice(0, nocc_*2));
      ctmp->copy_block(0, nocc_*2, coeff_->ndim(), nneg_/2, coeff_->slice(nocc_*2 + nvirtnr, nocc_*2 + nvirt_));
      ctmp->copy_block(0, nocc_*2 + nneg_/2, coeff_->ndim(), nneg_/2, coeff_->slice(nocc_*2 + nvirt_ + nvirtnr, nocc_*2 + nvirt_*2));
      *ctmp = *ctmp * *expa;
      auto ctmp2 = coeff_->copy();
      ctmp2->copy_block(0, 0, coeff_->ndim(), nocc_*2, ctmp->slice(0, nocc_*2));
      ctmp2->copy_block(0, nocc_*2 + nvirtnr, coeff_->ndim(), nneg_/2, ctmp->slice(nocc_*2, nocc_*2 +nneg_/2));
      ctmp2->copy_block(0, nocc_*2 + nvirtnr + nvirt_, coeff_->ndim(), nneg_/2, ctmp->slice(nocc_*2 + nneg_/2, ctmp->mdim()));
      coeff_ = make_shared<const ZMatrix>(*ctmp2);
    }
#endif
    // for next BFGS extrapolation
#ifdef BOTHSPACES
    *x *= *expa;
#else
    if (optimize_electrons) {
      *ele_x *= *expa;
    } else {
      *pos_x *= *expa;
    }
#endif

    // print energy
#ifdef BOTHSPACES
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());
#else
    if (optimize_electrons) {
      print_iteration(iter, 0, 0, ele_energy, gradient, timer.tick());
    } else {
      print_iteration(iter, 0, 0, pos_energy, gradient, timer.tick());
    }
#endif

#ifdef BOTHSPACES
    if (gradient < thresh_) break;
#else
    if (gradient < thresh_ && !optimize_electrons) pos_conv = true;
    if (gradient < thresh_ &&  optimize_electrons) ele_conv = true;
    optimize_electrons = optimize_electrons == true ? false : true;
    if (ele_conv && optimize_electrons) optimize_electrons = false;
    if (pos_conv && !optimize_electrons) optimize_electrons = true;
    if (pos_conv && ele_conv) break;
#endif
  }
//  if (ele_energy.size() > 0 && energy_.size() == 0) energy_.push_back(ele_energy.back());

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  if (nact_) {
    fci_->update(coeff_, /*restricted*/true);
    fci_->compute();
    fci_->compute_rdm12();
  }

}

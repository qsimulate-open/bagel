//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/asd_bfgs.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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

#include <src/scf/hf/fock.h>
#include <src/asd/orbital/asd_bfgs.h>
#include <src/asd/construct_asd.h>
#include <src/util/math/davidson.h>
#include <src/util/math/step_restrict_bfgs.h>
#include <src/util/math/hpw_diis.h>
//#include <src/util/math/bfgs.h>
#include <src/util/io/moldenout.h>

using namespace std;
using namespace bagel;

void ASD_BFGS::compute() {

  //BFGS parameters
  thresh_inter_ = idata_->get<double>("thresh_inter", 5.0e-4);
  thresh_intra_ = idata_->get<double>("thresh_intra", 5.0e-4);
  single_bfgs_ = idata_->get<bool>("single_bfgs", false);
  bool no_active = idata_->get<bool>("no_active", false);
  bool active_only = idata_->get<bool>("active_only", false);
  if (no_active || active_only) single_bfgs_ = false;

  vector<double> previous_energy(nstate_,0.0);

  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.
  shared_ptr<SRBFGS<ASD_RotFile>> bfgs;
//shared_ptr<BFGS<ASD_RotFile>> bfgs;
  shared_ptr<SRBFGS<ASD_RotFile>> intra_bfgs; //cloased-active, closed-virtual, active-virtual only
  shared_ptr<SRBFGS<ASD_RotFile>> inter_bfgs; //active-active only
  //shared_ptr<BFGS<ASD_RotFile>> intra_bfgs; //cloased-active, closed-virtual, active-virtual only
  //shared_ptr<BFGS<ASD_RotFile>> inter_bfgs; //active-active only
  pair<double,double> gradient_pair = make_pair(1.0,1.0);
  bool full = single_bfgs_ ? true : false;
  bool first_iteration = single_bfgs_ ? true : false;
  bool inter_first = single_bfgs_ ? false : true;
  bool intra_first = single_bfgs_ ? false : true;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  //allocate unitary rotation matrix
  auto x = make_shared<Matrix>(nbasis_, nbasis_);
  x->unit();
  auto intra_x = make_shared<Matrix>(nbasis_, nbasis_);
  intra_x->unit();
  auto inter_x = make_shared<Matrix>(nbasis_, nbasis_);
  inter_x->unit();
  vector<double> evals;
  vector<double> inter_evals;
  vector<double> intra_evals;

  cout << "     See asd_orbopt.log for further information on ASD output " << endl << endl;
  mute_stdcout();
  for (int iter = active_only ? 1 : 0; iter != max_iter_; ++iter) {

    bool inter = iter%2 == 0 ? true : false;

    //Perform ASD
    if (iter) dimer_->update_coeff(coeff_); //update coeff_ & integrals..
#if 1
    if (mpi__->rank() == 0) {
      string out_file = "orbital_" + to_string(iter) + ".molden";
      MoldenOut mfs(out_file);
      mfs << dimer_->sgeom();
      mfs << dimer_->sref();
    }
#endif
    Timer asd_time(0);
    auto asd = construct_ASD(idata_->get_child_optional("asd"), dimer_); //build CI-space with updated coeff
    asd_time.tick_print("ASD space construction");
    asd->compute();
    asd_time.tick_print("ASD");
    //get RDM
    rdm1_ = asd->rdm1_av();
    rdm2_ = asd->rdm2_av();
    asd_time.tick_print("RDM");
    //get energy
    energy_ = asd->energy();

    {
      // use state averaged energy to update trust radius
      double sa_en = 0.0;
      for (auto& i : asd->energy())
        sa_en += i;
      sa_en /= double((asd->energy()).size());
      if (full) evals.push_back(sa_en);
      else
        if (inter) inter_evals.push_back(sa_en);
        else       intra_evals.push_back(sa_en);
    }

    // compute one-body operators
    // * preparation
    const MatView ccoeff = coeff_->slice(0, nclosed_);

    // * core Fock operator
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<const Fock<1>>(geom_, hcore_, nullptr, ccoeff, /*store*/false, /*rhf*/true) : hcore_;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);

    // * active Fock operator
    // first make a weighted coefficient
    shared_ptr<Matrix> acoeff = coeff_->slice_copy(nclosed_, nocc_);
    shared_ptr<Matrix> rdm1_mat = rdm1_->rdm1_mat(/*nclose*/0);

    shared_ptr<Matrix> rdm1_scaled = rdm1_mat->copy();
    rdm1_scaled->sqrt();
    rdm1_scaled->delocalize();
    auto acoeffw = make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0)) * *rdm1_scaled); // such that C' * (1/2 D) C will be obtained.

    // then make a AO density matrix
    shared_ptr<const Matrix> afockao = make_shared<Fock<1>>(geom_, hcore_->clone(), nullptr, acoeffw, /*store*/false, /*rhf*/true);
    shared_ptr<const Matrix> afock = make_shared<Matrix>(*coeff_ % *afockao * *coeff_);

    // * Q_xr = 2(xs|tu)P_rs,tu (x=general, mo)
    auto qxr = Qvec(coeff_->mdim(), nact_, coeff_, nclosed_); //TODO: at the moment it uses internal rdm1_ and rdm2_

    //MCSCF Fock matrix (nact_,nact_) : active only
    shared_ptr<Matrix> mcfock = make_shared<Matrix>(*cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * *rdm1_mat
                                                    + *qxr->get_submatrix(nclosed_, 0, nact_, nact_) );

    //gradient
    auto grad = full ?  make_shared<ASD_RotFile>(nclosed_, nact_, nvirt_, rasA_, rasB_, true, true) :
                inter ? make_shared<ASD_RotFile>(nclosed_, nact_, nvirt_, rasA_, rasB_, true, false) :
                        make_shared<ASD_RotFile>(nclosed_, nact_, nvirt_, rasA_, rasB_, false, true);
    grad_vc(cfock, afock, grad); // grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
    grad_va(cfock, qxr, rdm1_mat, grad); // grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
    grad_ca(cfock, afock, qxr, rdm1_mat, grad); // grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
    grad_aa(mcfock, grad);

    // if this is the first time, set up the BFGS solver
    // BFGS and DIIS should start at the same time
    if (!single_bfgs_) {
      if (inter_first) { //if (iter == 0) {
        shared_ptr<const ASD_RotFile> inter_denom = compute_denom(cfock, afock, qxr, rdm1_mat, mcfock, true, false);
        inter_bfgs = make_shared<SRBFGS<ASD_RotFile>>(inter_denom);
//      inter_bfgs = make_shared<BFGS<ASD_RotFile>>(inter_denom);
        inter_first = false;
      }
      if (intra_first) { // (iter == 1) {
        cout << "intra_first" << endl;
        shared_ptr<const ASD_RotFile> intra_denom = compute_denom(cfock, afock, qxr, rdm1_mat, mcfock, false, true);
        intra_bfgs = make_shared<SRBFGS<ASD_RotFile>>(intra_denom);
//      intra_bfgs = make_shared<BFGS<ASD_RotFile>>(intra_denom);
        intra_first = false;
      }
    }

    if (full && first_iteration) {
      shared_ptr<const ASD_RotFile> denom = compute_denom(cfock, afock, qxr, rdm1_mat, mcfock, true, true);
      bfgs = make_shared<SRBFGS<ASD_RotFile>>(denom);
//    bfgs = make_shared<BFGS<ASD_RotFile>>(denom);
      first_iteration = false;
    }

    // extrapolation using BFGS
    Timer extrap(0);
    shared_ptr<ASD_RotFile> a;
    cout << " " << endl;
    cout << " -------  Step Restricted BFGS Extrapolation  ------- " << endl;
    if (full) {
      auto xcopy = x->log(8);
      auto xlog  = make_shared<ASD_RotFile>(xcopy, nclosed_, nact_, nvirt_, rasA_, rasB_, true, true);
      bfgs->check_step(evals, grad, xlog);
      a = bfgs->more_sorensen_extrapolate(grad, xlog);
//    a = bfgs->extrapolate(grad, xlog);
    } else {
      if (inter) {
        auto xcopy = inter_x->log(8);
        auto xlog  = make_shared<ASD_RotFile>(xcopy, nclosed_, nact_, nvirt_, rasA_, rasB_, true, false);
        inter_bfgs->check_step(inter_evals, grad, xlog); //, /*tight*/false, limited_memory);
        a = inter_bfgs->more_sorensen_extrapolate(grad, xlog);
//      a = inter_bfgs->extrapolate(grad, xlog);
      } else {
        auto xcopy = intra_x->log(8);
        auto xlog  = make_shared<ASD_RotFile>(xcopy, nclosed_, nact_, nvirt_, rasA_, rasB_, false, true);
        intra_bfgs->check_step(intra_evals, grad, xlog); //, /*tight*/false, limited_memory);
        a = intra_bfgs->more_sorensen_extrapolate(grad, xlog);
//      a = intra_bfgs->extrapolate(grad, xlog);
      }
    }
    cout << " ---------------------------------------------------- " << endl;
    extrap.tick_print("More-Sorensen/Hebden extrapolation");
    cout << " " << endl;

//  *a *= -1.0;
    a->print("Orbital rotation parameters");

    // restore the matrix from ASD_RotFile
    shared_ptr<const Matrix> amat = a->unpack<Matrix>();

    double max_rotation = 0.0;
    for (int i = 1; i != nbasis_; ++i)
      for (int j = 1; j != i; ++j)
        max_rotation = max(max_rotation, fabs(amat->element(j,i)));

    shared_ptr<Matrix> expa = amat->exp(100);
    expa->purify_unitary();

    // updating coefficients
    coeff_ = make_shared<const Coeff>(*coeff_ * *expa);

    // for next BFGS extrapolation
    if (full) *x *= *expa;
    else {
      if (inter) *inter_x *= *expa;
      else       *intra_x *= *expa;
    }
    // synchronization
    mpi__->broadcast(const_pointer_cast<Coeff>(coeff_)->data(), coeff_->size(), 0);

    // setting error of macro iteration
    const double gradient = grad->rms();

    //energy difference
    double delta_energy = 0.0;
    for (int i = 0; i != nstate_; ++i) {
      delta_energy = max(delta_energy, fabs(previous_energy[i] - energy_[i]));
    }
    copy(energy_.begin(), energy_.end(), previous_energy.begin() );

    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, max_rotation, delta_energy, timer.tick(), full ? 2 : static_cast<int>(inter));

    if (full && gradient < gradient_thresh_ && max_rotation < rotation_thresh_ && delta_energy < energy_thresh_) {
      rms_grad_ = gradient;
      cout << " " << endl;
      cout << "    * quasi-Newton optimization converged. *   " << endl << endl;
      mute_stdcout();
      break;
    }

    if (!full && inter) get<0>(gradient_pair) = gradient;
    else if (!full && !inter) get<1>(gradient_pair) = gradient;

    if (!full) {
      if (get<0>(gradient_pair) < thresh_inter_ && get<1>(gradient_pair) < thresh_intra_) {
        full = true;
        first_iteration = true;
        x->unit();
      }
    }
/*
    if (full && iter == 9) {
      if (gradient > 1.0e-4) {
        gradient_pair = make_pair(1.0,1.0);
        full = false;
        first_iteration = false;
        inter_first = true;
        intra_first = true;
        inter_x->unit();
        intra_x->unit();
      } else {
        first_iteration = true;
        x->unit();
      }
    }
*/

    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      if (rms_grad_ > gradient_thresh_) cout << "    * The calculation did NOT converge. *    " << endl;
      cout << "    * Max iteration reached during the quasi-Newton optimization. *     " << endl << endl;
    }
    mute_stdcout();
    if (active_only || no_active) ++iter;
  }
  resume_stdcout();
  // ============================
  // macro iteration to here
  // ============================
}

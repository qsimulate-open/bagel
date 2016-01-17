//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd/orbital/asd_bfgs.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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

#include <src/scf/hf/fock.h>
#include <src/asd/orbital/asd_bfgs.h>
#include <src/asd/construct_asd.h>
#include <src/util/math/step_restrict_bfgs.h>
#include <src/util/io/moldenout.h>

using namespace std;
using namespace bagel;

void ASD_BFGS::compute() {

  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  auto x = make_shared<Matrix>(nbasis_, nbasis_);
  x->unit();

  vector<double> evals;
  shared_ptr<SRBFGS<ASD_RotFile>> bfgs;
  double previous_energy = 0.0;
  bool macro = true;

  auto asd = construct_ASD(idata_->get_child_optional("asd"), dimer_, /*rdm=*/true);

  cout << endl << "     See asd_orbopt.log for further information on ASD output " << endl << endl;

  mute_stdcout();
  int iter = 0;
  int miter = -1;
  do {
    //Perform ASD
    Timer asd_time(0);
    if (iter) {
      dimer_->update_coeff(coeff_); //update coeff_ & integrals
      if (macro || miter < 0)
        asd = construct_ASD(idata_->get_child_optional("asd"), dimer_, /*rdm=*/true); //build CI-space with updated coeff
      else
        asd->update_dimer_and_fix_ci(dimer_); //fix ci coefficients
    }
    asd_time.tick_print("ASD space construction");
    asd->compute();
    asd_time.tick_print("ASD");

    //get RDM
    rdm1_ = asd->rdm1_av();
    rdm2_ = asd->rdm2_av();

    //get energy
    energy_ = asd->energy();

    //orbital
    if (print_orbital_ && mpi__->rank() == 0) {
      string out_file = "asd_orbital_iter_" + to_string(iter) + "_" + to_string(miter) + ".molden";
      MoldenOut mfs(out_file);
      mfs << dimer_->sgeom();
      mfs << dimer_->sref();
    }

    // use state averaged energy to update trust radius
    assert(asd->energy().size() == nstate_);
    double average_energy = 0.0;
    for (auto& i : asd->energy())
      average_energy += i;
    average_energy /= double((asd->energy()).size());

    evals.push_back(average_energy);

    // compute one-body Fock operators
    // preparation
    const MatView ccoeff = coeff_->slice(0, nclosed_);
    // core Fock operator
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<Fock<1>>(geom_, hcore_, nullptr, ccoeff, /*store*/false, /*rhf*/true) : hcore_;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    // active Fock operator
    // first make a weighted coefficient
    shared_ptr<const Matrix> acoeff = coeff_->slice_copy(nclosed_, nocc_);
    shared_ptr<const Matrix> rdm1_mat = rdm1_->rdm1_mat(/*nclose*/0);
    shared_ptr<Matrix> rdm1_scaled = rdm1_mat->copy();
    rdm1_scaled->sqrt();
    auto acoeffw = make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0)) * *rdm1_scaled); // such that C' * (1/2 D) C will be obtained.
    // then make a AO density matrix
    shared_ptr<const Matrix> afockao = make_shared<Fock<1>>(geom_, hcore_->clone(), nullptr, acoeffw, /*store*/false, /*rhf*/true);
    shared_ptr<const Matrix> afock = make_shared<Matrix>(*coeff_ % *afockao * *coeff_);
    // * Q_xr = 2(xs|tu)P_rs,tu (x=general, mo)
    auto qxr = Qvec(coeff_->mdim(), nact_, coeff_, nclosed_);
    //MCSCF Fock matrix (nact_,nact_) : active only
    shared_ptr<const Matrix> mcfock = make_shared<Matrix>((*cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * *rdm1_mat) + *qxr->get_submatrix(nclosed_, 0, nact_, nact_));

    //gradient evaluation
    auto grad = make_shared<ASD_RotFile>(nclosed_, nact_, nvirt_, rasA_, rasB_);
    grad_vc(cfock, afock, grad); // grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
    grad_va(cfock, qxr, rdm1_mat, grad); // grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
    grad_ca(cfock, afock, qxr, rdm1_mat, grad); // grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
    grad_aa(mcfock, grad);

    // if this is the first time (or first micro iteration), set up the BFGS solver
    if (miter < 0) {
      shared_ptr<const ASD_RotFile> denom = compute_denom(cfock, afock, qxr, rdm1_mat, mcfock);
      bfgs = make_shared<SRBFGS<ASD_RotFile>>(denom);
      miter = 0;
    }

    // extrapolation using BFGS
    Timer extrap(0);
    cout << " " << endl;
    cout << " -------  Step Restricted BFGS Extrapolation  ------- " << endl;
    auto xcopy = x->log(8);
    auto xlog  = make_shared<ASD_RotFile>(xcopy, nclosed_, nact_, nvirt_, rasA_, rasB_);
    bfgs->check_step(evals, grad, xlog);
    shared_ptr<const ASD_RotFile> a = bfgs->more_sorensen_extrapolate(grad, xlog);
    cout << " ---------------------------------------------------- " << endl;
    extrap.tick_print("More-Sorensen/Hebden extrapolation");
    cout << " " << endl;

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
    *x *= *expa;

    // synchronization
    mpi__->broadcast(const_pointer_cast<Coeff>(coeff_)->data(), coeff_->size(), 0);

    // setting error of macro iteration
    const double gradient = grad->rms();

    //energy difference
    const double delta_energy = average_energy - previous_energy;
    previous_energy = average_energy;

    resume_stdcout();
    print_iteration(iter, miter, 0, energy_, gradient, max_rotation, delta_energy, timer.tick());

    //user defined convergence criteria
    if (gradient < gradient_thresh_ && max_rotation < rotation_thresh_ && delta_energy < energy_thresh_ && miter == 0) {
      rms_grad_ = gradient;
      cout << " " << endl;
      cout << "    * quasi-Newton optimization converged. *   " << endl << endl;
      mute_stdcout();
      break;
    }

    if (macro) {
      if (miter == 0 && (gradient < fix_ci_begin_ || iter >= fix_ci_begin_iter_)) {
        // Start micro-iteration from next iteration
        iter++;
        miter = -1;
        x->unit();
        evals.clear();
        macro = false;
      } else
        iter++;
    } else {
      if (miter > 0 && gradient < fix_ci_thresh_) {
        if (fix_ci_finish_) {
          cout << endl << "    * quasi-Newton optimization converged with fixed monomer CI coefficients. *   " << endl << endl;
          break;
        }
        // Terminate micro-iteration, new monomer CI spaces will be formed
        iter++;
        miter = -1;
        x->unit();
        evals.clear();
      } else
        miter++;
    }

    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      cout << "    * The calculation did NOT converge. *    " << endl;
      cout << "    * Max iteration reached during the quasi-Newton optimization. *     " << endl << endl;
    }

    mute_stdcout();
  } while (iter < max_iter_-1);
  resume_stdcout();
  // ============================
  // macro iteration to here
  // ============================

  dimer_->update_coeff(coeff_);

  if (print_orbital_ && mpi__->rank() == 0) {
    string out_file = "asd_orbital_converged.molden";
    MoldenOut mfs(out_file);
    mfs << dimer_->sgeom();
    mfs << dimer_->sref();
  }

  if (semi_canonicalize_) {
    cout << "    * Form semi-canonical orbitals. *     " << endl << endl;
    coeff_ = semi_canonical_orb();
    dimer_->update_coeff(coeff_);
    if (print_orbital_ && mpi__->rank() == 0) {
      string out_file = "asd_orbital_semi_canonical.molden";
      MoldenOut mfs(out_file);
      mfs << dimer_->sgeom();
      mfs << dimer_->sref();
    }
  }

  // extra iteration for consistency
  asd = construct_ASD(idata_->get_child_optional("asd"), dimer_);
  asd->compute();

}

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: casbfgs1.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/multi/casscf/casbfgs.h>
#include <src/multi/casscf/qvec.h>
#include <src/scf/hf/fock.h>
#include <src/util/math/step_restrict_bfgs.h>
#include <src/prop/hyperfine.h>

using namespace std;
using namespace bagel;

void CASBFGS1::compute() {

  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.

  shared_ptr<SRBFGS<RotFile>> bfgs;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  auto x = make_shared<Matrix>(nbasis_, nbasis_);
  x->unit();
  shared_ptr<const Matrix> xstart;
  vector<double> evals;
  const int limited_memory = idata_->get<int>("limited_memory", 0);

  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {

    const shared_ptr<const Coeff> cold = coeff_;
    const shared_ptr<const Matrix> xold = x->copy();

    // first perform CASCI to obtain RDMs
    if (nact_) {
      if (iter) fci_->update(coeff_);
      Timer fci_time(0);
      fci_->compute();
      fci_->compute_rdm12();
      fci_time.tick_print("FCI and RDMs");
    }

    shared_ptr<Matrix> natorb_mat = x->clone();
    natorb_mat->unit();
    if (nact_) {
      // here make a natural orbitals and update coeff_. Closed and virtual orbitals remain canonical. Also, FCI::rdms are updated
      const pair<shared_ptr<Matrix>, VectorB> natorb = fci_->natorb_convert();
      coeff_ = update_coeff(coeff_, natorb.first);
      occup_ = natorb.second;
      if (natocc_) print_natocc();
      natorb_mat->copy_block(nclosed_, nclosed_, nact_, nact_, natorb.first);
    }

    auto sigma = make_shared<RotFile>(nclosed_, nact_, nvirt_);
    sigma->zero();

    // compute one-body operators
    Timer onebody;
    // * preparation
    const MatView ccoeff = coeff_->slice(0, nclosed_);
    // * core Fock operator
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<const Fock<1>>(geom_, hcore_, nullptr, ccoeff, /*store*/false, /*rhf*/true) : hcore_;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    // * active Fock operator
    shared_ptr<const Matrix> afock;
    if (nact_) {
      auto afockao = compute_active_fock(coeff_->slice(nclosed_, nocc_), fci_->rdm1_av());
      afock = make_shared<Matrix>(*coeff_ % *afockao * *coeff_);
    } else {
      afock = cfock->clone();
    }
    // * Q_xr = 2(xs|tu)P_rs,tu (x=general, mo)
    shared_ptr<const Qvec> qxr;
    if (nact_) {
      qxr = make_shared<const Qvec>(coeff_->mdim(), nact_, coeff_, nclosed_, fci_, fci_->rdm2_av());
    }

    // grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
    grad_vc(cfock, afock, sigma);
    // grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
    grad_va(cfock, qxr, sigma);
    // grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
    grad_ca(cfock, afock, qxr, sigma);

    // if this is the first time, set up the BFGS solver
    if (iter == 0) {
      // BFGS and DIIS should start at the same time
      shared_ptr<const RotFile> denom = compute_denom(cfock, afock, qxr);
      bfgs = make_shared<SRBFGS<RotFile>>(denom);
      const double trust_rad = idata_->get<double>("trust_radius", 0.4);
      bfgs->initiate_trust_radius(trust_rad);
    }
    onebody.tick_print("One body operators");

    // get energy
    if (nact_) {
      energy_ = fci_->energy();
      // use state averaged energy to update trust radius
      const double sa_en = blas::average(energy_);
      evals.push_back(sa_en);
    } else {
      const double en = (ccoeff % (*cfockao + *hcore_) * ccoeff).trace() + geom_->nuclear_repulsion();
      energy_ = {en};
      evals.push_back(en);
    }

    // extrapolation using BFGS
    Timer extrap(0);
    cout << " " << endl;
    cout << " -------  Step Restricted BFGS Extrapolation  ------- " << endl;
    *x *= *natorb_mat;
    auto xcopy = x->log(8);
    auto xlog  = make_shared<RotFile>(xcopy, nclosed_, nact_, nvirt_);
    bfgs->check_step(evals, sigma, xlog, /*tight*/false, limited_memory);
    shared_ptr<RotFile> a = bfgs->more_sorensen_extrapolate(sigma, xlog);
    cout << " ---------------------------------------------------- " << endl;
    extrap.tick_print("More-Sorensen/Hebden extrapolation");
    cout << " " << endl;

    // restore the matrix from RotFile
    shared_ptr<const Matrix> amat = a->unpack();
    shared_ptr<Matrix> expa = amat->exp(100);
    expa->purify_unitary();

    // updating coefficients
    coeff_ = make_shared<const Coeff>(*coeff_**expa);
    // for next BFGS extrapolation
    *x *= *expa;

    // synchronization
    mpi__->broadcast(const_pointer_cast<Coeff>(coeff_)->data(), coeff_->size(), 0);

    // setting error of macro iteration
    const double gradient = sigma->rms();

    resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    if (gradient < thresh_) {
      rms_grad_ = gradient;
      cout << " " << endl;
      cout << "    * quasi-Newton optimization converged. *   " << endl << endl;
      mute_stdcout();
      break;
    }

    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      if (rms_grad_ > thresh_) cout << "    * The calculation did NOT converge. *    " << endl;
      cout << "    * Max iteration reached during the quasi-Newton optimization. *     " << endl << endl;
      if (!conv_ignore_)
        throw runtime_error("CASSCF BFGS1 did not converge");
    }
    mute_stdcout();
  }
  resume_stdcout();
  // ============================
  // macro iteration to here
  // ============================

  // block diagonalize coeff_ in nclosed and nvirt
  coeff_ = semi_canonical_orb();

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  if (nact_) {
    fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
  }

  // calculate the HFCCs
  if (do_hyperfine_ && !geom_->external() && nstate_ == 1) {
    HyperFine hfcc(geom_, spin_density(), fci_->det()->nspin(), "CASSCF");
    hfcc.compute();
  }
}

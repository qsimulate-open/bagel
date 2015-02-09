//
// BAGEL - Parallel electron correlation program.
// Filename: multi/rasscf/bfgs.cc
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


#include <src/multi/rasscf/bfgs.h>
#include <src/multi/rasscf/qvec.h>
#include <src/scf/hf/fock.h>
#include <src/util/math/davidson.h>
#include <src/util/math/step_restrict_bfgs.h>
#include <src/util/math/hpw_diis.h>

using namespace std;
using namespace bagel;

void RASBFGS::compute() {

  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.

  shared_ptr<SRBFGS<RASRotFile>> bfgs;

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

    // first perform RASCI to obtain RDMs
    if (nact_) {
      if (iter) rasci_->update(coeff_);
      Timer rasci_time(0);
      rasci_->compute();
      rasci_->compute_rdm12();
      rasci_time.tick_print("RASCI and RDMs");
    }
 
    cout << "RASBFGS:: check 1" << endl;

//TODO: temporary
    const MatView cdata = coeff_->slice(nclosed_, nclosed_+nact_);
    half_ = geom_->df()->compute_half_transform(cdata);

    cout << "RASBFGS:: check 2" << endl;
    auto sigma = make_shared<RASRotFile>(nclosed_, nact_, nvirt_, ras_);
    sigma->zero();

    cout << "RASBFGS:: check 3" << endl;
    // compute one-body operators
    Timer onebody;
    // * preparation
    const MatView ccoeff = coeff_->slice(0, nclosed_);
    // * core Fock operator
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<const Fock<1>>(geom_, hcore_, nullptr, ccoeff, /*store*/false, /*rhf*/true) : hcore_;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    cout << "RASBFGS:: check 4" << endl;
    // * active Fock operator
    // first make a weighted coefficient
//  shared_ptr<Matrix> acoeff;
//  if (nact_) {
//    acoeff = coeff_->slice_copy(nclosed_, nocc_);
//    for (int i = 0; i != nact_; ++i)
//      blas::scale_n(sqrt(occup_[i]/2.0), acoeff->element_ptr(0, i), acoeff->ndim());
//  }
    shared_ptr<Matrix> acoeff = coeff_->slice_copy(nclosed_, nocc_);
    shared_ptr<Matrix> rdm1_mat = rasci_->rdm1_av()->rdm1_mat(/*nclose*/0);
    cout << "RASBFGS:: check 5" << endl;

    shared_ptr<Matrix> rdm1_scaled = rdm1_mat->copy();
    rdm1_scaled->sqrt();
    rdm1_scaled->delocalize();
    auto acoeffw = make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0)) * *rdm1_scaled);
    cout << "RASBFGS:: check 6" << endl;

    // then make a AO density matrix
    shared_ptr<const Matrix> afock;
    if (nact_) {
      auto afockao = make_shared<Fock<1>>(geom_, hcore_, nullptr, acoeffw, /*store*/false, /*rhf*/true);
      afock = make_shared<Matrix>(*coeff_ % (*afockao - *hcore_) * *coeff_);
    } else {
      afock = cfock->clone();
    }
    cout << "RASBFGS:: check 7" << endl;
    // * Q_xr = 2(xs|tu)P_rs,tu (x=general, mo)
    shared_ptr<const RASQvec> qxr;
    if (nact_) {
      qxr = make_shared<const RASQvec>(coeff_->mdim(), nact_, coeff_, nclosed_, rasci_, rasci_->rdm2_av());
    }
    cout << "RASBFGS:: check 8" << endl;

    //mcfock(nact_,nact_)
    shared_ptr<Matrix> mcfock = make_shared<Matrix>(*cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * *rdm1_mat
                                                    + *qxr->get_submatrix(nclosed_, 0, nact_, nact_) );
    cout << "RASBFGS:: check 9" << endl;

    // grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
    grad_vc(cfock, afock, sigma);
    // grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
    grad_va(cfock, qxr, rdm1_mat, sigma);
    // grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
    grad_ca(cfock, afock, qxr, rdm1_mat, sigma);
    cout << "RASBFGS:: check 10" << endl;

    grad_aa12(mcfock, sigma);
    grad_aa13(mcfock, sigma);
    grad_aa23(mcfock, sigma);
    cout << "RASBFGS:: check 11" << endl;

    // if this is the first time, set up the BFGS solver
    if (iter == 0) {
      // BFGS and DIIS should start at the same time
      shared_ptr<const RASRotFile> denom = compute_denom(cfock, afock, qxr, rdm1_mat, mcfock);
      bfgs = make_shared<SRBFGS<RASRotFile>>(denom);
    }
    onebody.tick_print("One body operators");

    cout << "RASBFGS:: check 12" << endl;

    // get energy
    if (nact_) {
      energy_ = rasci_->energy();
      // use state averaged energy to update trust radius
      const double sa_en = blas::average(energy_);
      evals.push_back(sa_en);
    } else {
      const double en = (ccoeff % (*cfockao + *hcore_) * ccoeff).trace() + geom_->nuclear_repulsion();
      energy_ = {en};
      evals.push_back(en);
    }

    cout << "RASBFGS:: check 13" << endl;
    // extrapolation using BFGS
    Timer extrap(0);
    cout << " " << endl;
    cout << " -------  Step Restricted BFGS Extrapolation  ------- " << endl;
//  *x *= *natorb_mat;
    auto xcopy = x->log(8);
    auto xlog  = make_shared<RASRotFile>(xcopy, nclosed_, nact_, nvirt_, ras_);
    bfgs->check_step(evals, sigma, xlog, /*tight*/false, limited_memory);
    shared_ptr<RASRotFile> a = bfgs->more_sorensen_extrapolate(sigma, xlog);
    cout << " ---------------------------------------------------- " << endl;
    extrap.tick_print("More-Sorensen/Hebden extrapolation");
    cout << " " << endl;

    // restore the matrix from RASRotFile
    shared_ptr<const Matrix> amat = a->unpack<Matrix>();
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
    }
    mute_stdcout();
  }
  resume_stdcout();
  assert(false);
  // ============================
  // macro iteration to here
  // ============================

  // block diagonalize coeff_ in nclosed and nvirt
//coeff_ = semi_canonical_orb();

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  if (nact_) {
    rasci_->update(coeff_);
    rasci_->compute();
    rasci_->compute_rdm12();
  }
}


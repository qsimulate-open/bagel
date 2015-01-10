//
// BAGEL - Parallel electron correlation program.
// Filename: casbfgs.cc
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


#include <src/asdscf/bfgs2.h> //TODO bfgs2->bfgs
#include <src/math/davidson.h>
#include <src/math/step_restrict_bfgs.h>
#include <src/math/hpw_diis.h>

#include <src/asd/construct_asd.h>

using namespace std;
using namespace bagel;

void ASDBFGS2::compute() {

  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.

  shared_ptr<SRBFGS<ASDRotFile2>> bfgs;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  //allocate unitary rotation matrix
  auto x = make_shared<Matrix>(nbasis_, nbasis_);
  x->unit();
  shared_ptr<const Matrix> xstart;
  vector<double> evals;

  auto asd = construct_ASD(asdinput_, dimer_);
  rdm1_ = make_shared<RDM<1>>(nact_);
  rdm2_ = make_shared<RDM<2>>(nact_);

//mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {

    const shared_ptr<const Coeff> cold = coeff_;
    const shared_ptr<const Matrix> xold = x->copy();

    // first perform CASCI to obtain RDMs
//  if (iter) fci_->update(coeff_);
//  fci_->compute();
//  fci_->compute_rdm12();
//  // get energy
//  energy_ = fci_->energy();

    //Perform ASD
    if (iter) {
      //update coeff_ & integrals..
      cout << "BFGS: update coeff" << endl;
      coeff_->print();
      dimer_->update_coeff(coeff_);
      //build CI-space with updated coeff
      asd = construct_ASD(asdinput_, dimer_);
    }
    cout << "BFGS: ASD.." << endl;
    asd->compute();
    cout << "BFGS: ASD done.." << endl;
    //get RDM
    rdm1_ = asd->rdm1();
    rdm2_ = asd->rdm2();
    //get energy
    energy_ = asd->energy();

  //{
  //  // use state averaged energy to update trust radius
  //  double sa_en = 0.0;
  //  for (auto& i : fci_->energy())
  //    sa_en += i;
  //  sa_en /= double((fci_->energy()).size());
  //  evals.push_back(sa_en);
  //}
//ONLY GROUND STATE
    evals.push_back(energy_[0]);
    cout << "BFGS: evals done.." << endl;
    
//Disable Natural orbital
/*
    shared_ptr<Matrix> natorb_mat = x->clone();
    {
      // here make a natural orbitals and update coeff_. Closed and virtual orbitals remain canonical. Also, FCI::rdms are updated
      shared_ptr<const Matrix> natorb = form_natural_orbs();
      natorb_mat->unit();
      natorb_mat->copy_block(nclosed_, nclosed_, nact_, nact_, natorb);
    }
    cout << "BFGS: natural orbital done.." << endl;
*/
//but get half_
    const MatView cdata = coeff_->slice(nclosed_, nclosed_+nact_);
    half_ = geom_->df()->compute_half_transform(cdata);
//END

    auto sigma = make_shared<ASDRotFile2>(nclosed_, nact_, nvirt_, nactA_, nactB_ );
    sigma->zero();

    // compute one-body operators
    // * preparation
    const MatView ccoeff = coeff_->slice(0, nclosed_);
    cout << "BFGS: 1.." << endl;



    // * core Fock operator
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<const Fock<1>>(geom_, hcore_, nullptr, ccoeff, /*store*/false, /*rhf*/true) : hcore_;
    cout << "BFGS: 2.." << endl;
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    cout << "BFGS: Fock(closed) done..: size = " << cfock->ndim() << " x " << cfock->mdim() << endl;


    // * active Fock operator
    // first make a weighted coefficient
    shared_ptr<Matrix> acoeff = coeff_->slice_copy(nclosed_, nocc_);
    shared_ptr<Matrix> rdm1_mat = rdm1_->rdm1_mat(/*nclose*/0);

    shared_ptr<Matrix> rdm1_scaled = rdm1_mat->copy();
    rdm1_scaled->sqrt();
    rdm1_scaled->delocalize();
    auto acoeffw = make_shared<Matrix>(*acoeff * (1.0/sqrt(2.0)) * *rdm1_scaled); // such that C' * (1/2 D) C will be obtained.

//  for (int i = 0; i != nact_; ++i)
//    blas::scale_n(sqrt(occup_[i]/2.0), acoeff->element_ptr(0, i), acoeff->ndim());
    // then make a AO density matrix
    shared_ptr<const Matrix> afockao = make_shared<Fock<1>>(geom_, hcore_->clone(), nullptr, acoeffw, /*store*/false, /*rhf*/true);
    shared_ptr<const Matrix> afock = make_shared<Matrix>(*coeff_ % *afockao * *coeff_);
    cout << "BFGS: Fock(active) done..: size = " << afock->ndim() << " x " << afock->mdim() << endl;



    // * Q_xr = 2(xs|tu)P_rs,tu (x=general, mo)
//  auto qxr = make_shared<const Qvec>(coeff_->mdim(), nact_, coeff_, nclosed_, fci_, fci_->rdm2_av());
//ADDED
    auto qxr = Qvec(coeff_->mdim(), nact_, coeff_, nclosed_); //uses internal rdm1_ and rdm2_
    cout << "BFGS: Qvec done.." << endl;


    //mcfock
    shared_ptr<Matrix> mcfock = make_shared<Matrix>(*cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * *rdm1_mat
                                                    + *qxr->get_submatrix(nclosed_, 0, nact_, nact_) );
    cout << "MC Fock Matrix: symmetric check:" << check_symmetric(mcfock) << endl;
    {
      shared_ptr<Matrix> p = make_shared<Matrix>(*cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * *rdm1_mat );
      shared_ptr<Matrix> q = make_shared<Matrix>( *qxr->get_submatrix(nclosed_, 0, nact_, nact_) );
      cout << "TEST P Matrix: symmetric check:" << check_symmetric(p) << endl;
      cout << "TEST Q Matrix: symmetric check:" << check_symmetric(q) << endl;
    }
  


    // grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
    grad_vc(cfock, afock, sigma);
    cout << "BFGS: Gradient vc done.." << endl;
    // grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
    grad_va(cfock, qxr, rdm1_mat, sigma);
    cout << "BFGS: Gradient va done.." << endl;
    // grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
    grad_ca(cfock, afock, qxr, rdm1_mat, sigma);
    cout << "BFGS: Gradient ca done.." << endl;
    grad_aa(mcfock, sigma);
    cout << "BFGS: Gradient aa done.." << endl;

    // if this is the first time, set up the BFGS solver
    if (iter == 0) {
      // BFGS and DIIS should start at the same time
      shared_ptr<const ASDRotFile2> denom = compute_denom(cfock, afock, qxr, rdm1_mat, mcfock);
      bfgs = make_shared<SRBFGS<ASDRotFile2>>(denom);
    }

    // extrapolation using BFGS
    cout << " " << endl;
    cout << " -------  Step Restricted BFGS Extrapolation  ------- " << endl;
//  *x *= *natorb_mat;
    auto xcopy = x->log(8);
    auto xlog  = make_shared<ASDRotFile2>(xcopy, nclosed_, nact_, nvirt_, nactA_, nactB_ );
    bfgs->check_step(evals, sigma, xlog);
    shared_ptr<ASDRotFile2> a = bfgs->more_sorensen_extrapolate(sigma, xlog);
    cout << " ---------------------------------------------------- " << endl << endl;

    // restore the matrix from ASDRotFile2
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

//  resume_stdcout();
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    if (gradient < thresh_) {
      rms_grad_ = gradient;
      cout << " " << endl;
      cout << "    * quasi-Newton optimization converged. *   " << endl << endl;
//    mute_stdcout();
      break;
    }

    if (iter == max_iter_-1) {
      rms_grad_ = gradient;
      cout << " " << endl;
      if (rms_grad_ > thresh_) cout << "    * The calculation did NOT converge. *    " << endl;
      cout << "    * Max iteration reached during the quasi-Newton optimization. *     " << endl << endl;
    }
//  mute_stdcout();
  }
  cout << "BFGS macro end" << endl;
//resume_stdcout();
  // ============================
  // macro iteration to here
  // ============================

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
//SKIPPED
//fci_->update(coeff_);
//fci_->compute();
//fci_->compute_rdm12();
}


shared_ptr<const ASDRotFile2> ASDBFGS2::compute_denom(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> mcfock) const {
  auto out = make_shared<ASDRotFile2>(nclosed_, nact_, nvirt_, nactA_, nactB_ );
//const double tiny = 1.0e-15;

  shared_ptr<Matrix> cfockd;
  if (nact_) {
    cfockd = make_shared<Matrix>(*cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * *rdm1);
    //TODO check symmetric??
  }

  // ia part (4.7a)
  if (nvirt_ && nclosed_) {
    double* target = out->ptr_vc();
    for (int i = 0; i != nclosed_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        *target++ = 4.0*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_)) - 4.0*(cfock->element(i,i)+afock->element(i,i));
      }
    }
  }
  // ra part (4.7b)
  if (nvirt_ && nact_) {
    double* target = out->ptr_va();
    for (int i = 0; i != nact_; ++i) {
    //if (occup_[i] < tiny) continue;
      for (int j = 0; j != nvirt_; ++j) {
    //  *target++ = 2.0*occup_[i]*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_))
    //            - 2.0*occup_[i]*cfock->element(i+nclosed_, i+nclosed_) - 2.0*qxr->element(i+nclosed_, i);
        *target++ = 2.0*rdm1->element(i,i)*(cfock->element(j+nocc_, j+nocc_)+afock->element(j+nocc_, j+nocc_))
                  - 2.0*cfockd->element(i,i) - 2.0*qxr->element(i+nclosed_, i);
      }
    }
  }
  // it part (4.7c)
  if (nclosed_ && nact_) {
    double* target = out->ptr_ca();
    for (int i = 0; i != nact_; ++i) {
//    if (occup_[i] < tiny) continue;
      for (int j = 0; j != nclosed_; ++j) {
//      *target++ = 4.0*(cfock->element(i+nclosed_, i+nclosed_)+afock->element(i+nclosed_, i+nclosed_) - cfock->element(j,j) - afock->element(j,j))
//                + 2.0*occup_[i]*(cfock->element(j,j)+afock->element(j,j)) - 2.0*occup_[i]*cfock->element(i+nclosed_, i+nclosed_) - 2.0*qxr->element(i+nclosed_, i);
        *target++ = 4.0*(cfock->element(i+nclosed_, i+nclosed_)+afock->element(i+nclosed_, i+nclosed_) - cfock->element(j,j) - afock->element(j,j))
                  + 2.0*rdm1->element(i,i)*(cfock->element(j,j)+afock->element(j,j)) - 2.0*cfockd->element(i,i) - 2.0*qxr->element(i+nclosed_, i);
      }
    }
  }
  // tu part
  if (nact_) {
    double* target = out->ptr_aa();
    for (int i = nactA_; i != nact_; ++i) { //B
      for (int j = 0; j != nactA_; ++j) { //A
        *target++ = 2.0*( 
                    - mcfock->element(j,j) - mcfock->element(i,i) 
                    + rdm1->element(j,j)*cfock->element(i,i) + rdm1->element(i,i)*cfock->element(j,j) 
                    - rdm1->element(i,j)*cfock->element(j,i) - rdm1->element(j,i)*cfock->element(i,j)
                    + 2.0*afock->element(i,i) + 2.0*afock->element(j,j)
                    );
      }
    }
  }

  const double thresh = 1.0e-8;
  for (int i = 0; i != out->size(); ++i)
    if (fabs(out->data(i)) < thresh) {
      out->data(i) = 1.0e10;
    }
  return out;
}


// grad(a/i) (eq.4.3a): 4(cfock_ai+afock_ai)
void ASDBFGS2::grad_vc(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<ASDRotFile2> sigma) const {
  if (!nvirt_ || !nclosed_) return;
  double* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_; ++i, target += nvirt_) {
    daxpy_(nvirt_, 4.0, cfock->element_ptr(nocc_,i), 1, target, 1);
    daxpy_(nvirt_, 4.0, afock->element_ptr(nocc_,i), 1, target, 1);
  }
}


// grad(a/t) (eq.4.3b): 2cfock_au gamma_ut + q_at
void ASDBFGS2::grad_va(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> qxr, shared_ptr<Matrix> rdm1, shared_ptr<ASDRotFile2> sigma) const {
  if (!nvirt_ || !nact_) return;
  dgemm_("N", "T", nvirt_, nact_, nact_, 2.0, cfock->element_ptr(nocc_,nclosed_), cfock->ndim(), rdm1->data(), rdm1->ndim(), 0.0, sigma->ptr_va(), nvirt_);
  double* target = sigma->ptr_va();
  for (int i = 0; i != nact_; ++i, target += nvirt_) {
//  daxpy_(nvirt_, 2.0*occup_[i], cfock->element_ptr(nocc_, i+nclosed_), 1, target, 1);
    daxpy_(nvirt_, 2.0, qxr->element_ptr(nocc_, i), 1, target, 1);
  }
}


// grad(r/i) (eq.4.3c): 4(cfock_ri+afock_ri) - 2cfock_iu gamma_ur - qxr_ir
void ASDBFGS2::grad_ca(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr, shared_ptr<Matrix> rdm1, shared_ptr<ASDRotFile2> sigma) const {
  if (!nclosed_ || !nact_) return;
  {
    double* target = sigma->ptr_ca();
    for (int i = 0; i != nact_; ++i, target += nclosed_) {
    //daxpy_(nclosed_, 4.0-2.0*occup_[i], cfock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, 4.0, cfock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, 4.0, afock->element_ptr(0,nclosed_+i), 1, target, 1);
      daxpy_(nclosed_, -2.0, qxr->element_ptr(0, i), 1, target, 1);
    }
    //-2 cfock_iu * D_ur
    dgemm_("T", "N", nclosed_, nact_, nact_, -2.0, cfock->element_ptr(nclosed_,0), cfock->ndim(), rdm1->data(), rdm1->ndim(), 1.0, sigma->ptr_ca(), nclosed_);
  }
}

// grad(t/t)
void ASDBFGS2::grad_aa(shared_ptr<const Matrix> mcfock, shared_ptr<ASDRotFile2> sigma) const {
  if (!nact_) return;
  double* target = sigma->ptr_aa();
  for (int i = 0; i != nactB_; ++i) { //B
    for (int j = 0; j != nactA_; ++j, ++target) { //A
      *target = 2.0*(mcfock->element(j,i+nactA_) - mcfock->element(i+nactA_,j));
    //*target = - mcfock->element(j,i+nactA_) + mcfock->element(i+nactA_,j);
    }
  }
}

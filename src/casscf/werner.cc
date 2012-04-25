//
// Newint - Parallel electron correlation program.
// Filename: werner.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//



#include <src/casscf/werner.h>
#include <src/casscf/jvec.h>
#include <src/scf/matrix1e.h>
#include <src/casscf/rotfile.h>
#include <src/util/aughess.h>
#include <src/util/linear.h>

using namespace std;

static const double cps = static_cast<double>(CLOCKS_PER_SEC);

void WernerKnowles::compute() {
  const string indent = "  ";

  cout << indent << "=== CASSCF iteration (" + geom_->basisfile() + ") ===" << endl << endl;

  // initializing Hcore matrix (redundant copy, but I can live with it).
  {
    shared_ptr<Hcore> hc(new Hcore(geom_));
    shared_ptr<Fock<1> > fc(new Fock<1>(geom_, hc));
    hcore_ = fc;
  }

  // All the equation numbers refer to those in Werner and Knowles, J. Chem. Phys. 82, 5053 (1985).
  // macro iteration
  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {

    int start = ::clock();

    // performs FCI, which also computes one-index transformed integrals
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    vector<double> energy = fci_->energy();

    // from natural orbitals from FCI::rdm1_av_ and set appropriate quantities.
#if 0
    form_natural_orbs();
#else
    occup_ = fci_->rdm1_av()->diag();
#endif

    shared_ptr<Jvec> jvec(new Jvec(fci_, ref_->coeff(), nclosed_, nact_, nvirt_));

    // denominator
    shared_ptr<Matrix1e> denom;
    {
      shared_ptr<Matrix1e> f;
      shared_ptr<QFile>    fact, factp, gaa;
      shared_ptr<RotFile> denom_;
      one_body_operators(f, fact, factp, gaa, denom_, false);
      denom = denom_->unpack(geom_, 1.0e-1);
    }

    // start with U=1
    shared_ptr<Matrix1e> U(new Matrix1e(geom_));
    {
      U->unit();
      shared_ptr<Matrix1e> bvec = compute_bvec(jvec, U, ref_->coeff());
      shared_ptr<Matrix1e> residual(new Matrix1e(*bvec - *bvec->transpose()));
      for (int i = 0; i != residual->size(); ++i) residual->data(i) /=  max(std::abs(denom->data(i)),0.1);
      U = residual->exp();
      U->purify_unitary();
    }


    shared_ptr<Matrix1e> bvec0;

    int miter;
    double error;
    for (miter = 0; miter != max_micro_iter_; ++miter) {

      // compute initial B (Eq. 19)
      shared_ptr<Matrix1e> bvec = compute_bvec(jvec, U, ref_->coeff());
      if (miter == 0) bvec0 = bvec;

      { // for debug
        shared_ptr<Matrix1e> t(new Matrix1e(*U));
        for (int i = 0; i != t->ndim(); ++i) t->element(i,i) -= 1.0;
        const double energy_micro = (bvec0->ddot(*t)*0.5+ bvec->ddot(*t)*0.5)+ energy[0];
        cout << "micro energy " << setw(20) << setprecision(10) << energy_micro << " " << bvec0->ddot(*t) << " " << bvec->ddot(*t) << endl;
      }

      // compute gradient
      shared_ptr<Matrix1e> grad(new Matrix1e(*U%*bvec-*bvec%*U));
grad->print("grad", 18);
      grad->purify_redrotation(nclosed_,nact_,nvirt_);
      const double error_micro = grad->ddot(*grad)/grad->size();
      if (miter == 0) error = error_micro;
      cout << "   -- iter " << setw(4) << miter << "  residual norm : " << setw(20) << setprecision(10) << error_micro << endl;

      if (error_micro < thresh_micro_) break;

      // initializing a Davidson manager
#if 1
      AugHess<Matrix1e> solver(max_mmicro_iter_+1, grad);
#else
      Linear<Matrix1e> solver(max_mmicro_iter_, grad);
#endif

      // update C = 1/2(A+A^dagger) = 1/2(U^dagger B + B^dagger U)
      shared_ptr<Matrix1e> C(new Matrix1e((*U % *bvec + *bvec % *U)*0.5));

      // initial dR value.
      shared_ptr<Matrix1e> dR(new Matrix1e(*grad));
      for (int i = 0; i != dR->size(); ++i) dR->data(i) /=  max(std::abs(denom->data(i)),0.1);

      // solving Eq. 30 of Werner and Knowles
      // fixed U, solve for dR.
      for (int mmiter = 0; mmiter != max_mmicro_iter_; ++mmiter) {

        const int mstart = ::clock();
        // compute U dR
        shared_ptr<Matrix1e> UdR(new Matrix1e(*U**dR));

        // update B
        shared_ptr<Matrix1e> new_bvec = compute_bvec(jvec, UdR, UdR, ref_->coeff());

        // compute  C dR + dR C
        shared_ptr<Matrix1e> dRA(new Matrix1e(*C**dR+*dR**C));

        // compute U^dagger B - B^dagger U
        // compute Eq.29
        shared_ptr<Matrix1e> sigma(new Matrix1e(*U%*new_bvec-*new_bvec%*U-*dRA));
        sigma->purify_redrotation(nclosed_,nact_,nvirt_);

        shared_ptr<Matrix1e> residual = solver.compute_residual(dR, sigma);

        const double error_mmicro = residual->ddot(*residual) / residual->size();
        const double mic_energy = 0.0; // TODO

        const int mend = ::clock();
#if 1
        if (mmiter == 0) cout << endl << "     == micro iteration == " << endl;
        cout << setw(10) << mmiter << "   " << setw(20) << setprecision(12)
             << setw(20) << setprecision(14) << error_mmicro << " " << setprecision(2) << (mend - mstart)/cps << endl;
#endif

        if (error_mmicro < thresh_mmicro_) { cout << endl; break; }
        if (mmiter+1 == max_mmicro_iter_) { cout << endl; break; }

        // update dR;
        for (int i = 0; i != residual->size(); ++i) residual->data(i) /=  max(std::abs(denom->data(i)),0.1);
        solver.orthog(residual);
        dR = residual;
        dR->purify_redrotation(nclosed_,nact_,nvirt_);

      }
      dR = solver.civec();
      shared_ptr<Matrix1e> dU = dR->exp();
      dU->purify_unitary();
      *U *= *dU;

    }

#if 1
((*U-*(U->transpose()))*0.5).print("uu", 18);
shared_ptr<Matrix1e> bvec2 = compute_bvec(jvec, U, ref_->coeff());
(*U%*bvec2-*bvec2%*U).print("residual",18);
#endif

    shared_ptr<Coeff> newcc(new Coeff(*ref_->coeff() * *U));
    ref_->set_coeff(newcc);

    int end = ::clock();
    if (nstate_ != 1 && iter) cout << endl;
    for (int i = 0; i != nstate_; ++i) {
      resume_stdcout();
      cout << indent << setw(5) << iter << setw(3) << i << setw(4) << miter
                                        << setw(15) << fixed << setprecision(8) << energy[i] << "   "
                                        << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0) << fixed << setw(10) << setprecision(2)
                                        << (end - start)/cps << endl;
      mute_stdcout();
    }
    if (error < thresh_) break;

  }
  resume_stdcout();

}



// compute B according to Eq. (19).
// B = 2[h_rs U_sj D_ji + (rs|D)Usj <D|ji> + 2(rk|D)(D|ls)T_sj D_jl,ik]
shared_ptr<Matrix1e> WernerKnowles::compute_bvec(shared_ptr<Jvec> jvec, shared_ptr<Matrix1e> u, const shared_ptr<Coeff> cc) {
  shared_ptr<Matrix1e> t(new Matrix1e(*u));
  for (int i = 0; i != u->ndim(); ++i) t->element(i,i) -= 1.0;
  return compute_bvec(jvec, u, t, cc);
}


shared_ptr<Matrix1e> WernerKnowles::compute_bvec(shared_ptr<Jvec> jvec,
                                                 shared_ptr<Matrix1e> u, shared_ptr<Matrix1e> t, const shared_ptr<Coeff> cc) {
  shared_ptr<DensityFit> df = geom_->df();
  const int naux = df->naux();
  const int nbas = df->nbasis();

// TODO make sure if this works
if (nbasis_ != nbas) throw runtime_error("I should examine this case...");

  shared_ptr<Matrix1e> out(new Matrix1e(geom_));
  {
    shared_ptr<Matrix1e> hcore_mo_(new Matrix1e(*cc % *hcore_ * *cc));

    // first term
    Matrix1e all1(geom_);
    for (int i = 0; i != nclosed_; ++i) all1.element(i,i) = 2.0;
    for (int i = 0; i != nact_; ++i) dcopy_(nact_, fci_->rdm1_av()->data()+nact_*i, 1, all1.element_ptr(nclosed_, i+nclosed_), 1);
    shared_ptr<Matrix1e> buf(new Matrix1e(geom_));
    dgemm_("N", "N", nbasis_, nocc_, nbasis_, 1.0, hcore_mo_->data(), nbasis_, u->data(), nbas, 0.0, buf->data(), nbas);
    dgemm_("N", "N", nbasis_, nocc_, nocc_, 2.0, buf->data(), nbas, all1.data(), nbas, 0.0, out->data(), nbas);
  }

  unique_ptr<double[]> tmp(new double[nbas*nbas]);
  // second term
  {
    Matrix1e Umat(*cc * *u);
    shared_ptr<DF_Half> half = df->compute_half_transform(Umat.data(), nocc_);
    half->form_2index(tmp, jvec->jvec(), 2.0, 0.0);
  }

  // third term
  if (t->norm() > 1.0e-15) {
    Matrix1e Tmat(*cc * *t);
    shared_ptr<DF_Full> full = jvec->half()->compute_second_transform(Tmat.data(), nocc_)->apply_2rdm(jvec->rdm2_all());
    jvec->half()->form_2index(tmp, full, 4.0, 1.0);
  }

  dgemm_("T", "N", nbasis_, nocc_, nbas, 1.0, cc->data(), nbas, tmp.get(), nbas, 1.0, out->data(), nbasis_);

  return out;
}


//
// BAGEL - Parallel electron correlation program.
// Filename: werner.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <src/casscf/werner.h>
#include <src/casscf/jvec.h>
#include <src/scf/matrix1e.h>
#include <src/casscf/rotfile.h>
#include <src/util/aughess.h>
#include <src/util/linear.h>
#include <src/casscf/jkop.h>
#include <src/util/bfgs.h>

using namespace std;
using namespace bagel;

static const vector<double> zero(1,0.0);

//#define SOLVER Linear
#define SOLVER AugHess

void WernerKnowles::compute() {

  // All the equation numbers refer to those in Werner and Knowles, J. Chem. Phys. 82, 5053 (1985).
  // macro iteration
  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer timer;

    // performs FCI, which also computes one-index transformed integrals
    fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    energy_ = fci_->energy();

    // from natural orbitals from FCI::rdm1_av_ and set appropriate quantities.
    form_natural_orbs();
    shared_ptr<Jvec> jvec(new Jvec(fci_, coeff_, nclosed_, nact_, nvirt_));

    // denominator

    // start with U=1
    shared_ptr<Matrix> U(new Matrix(geom_->nbasis(), geom_->nbasis()));
    U->unit();


    int miter, tcount = 0;
    double error = 0.0;
    for (miter = 0; miter != max_micro_iter_; ++miter) {

      // compute initial B (Eq. 19)
      shared_ptr<Matrix> bvec = compute_bvec(jvec, U, coeff_);

      // compute gradient
      shared_ptr<Matrix> grad(new Matrix(*U%*bvec-*bvec%*U));
      grad->purify_redrotation(nclosed_,nact_,nvirt_);

      const double error_micro = grad->ddot(*grad)/grad->size();
      if (miter == 0) error = error_micro;
      cout << "   -- iter " << setw(4) << miter << "  residual norm : " << setw(20) << setprecision(10) << error_micro << endl;

      if (error_micro < thresh_micro_) break;

      // initializing a linear solver (SOLVER is defined above in this file)
      SOLVER<Matrix> solver(max_mmicro_iter_+1, grad);

      // update C = 1/2(A+A^dagger) = 1/2(U^dagger B + B^dagger U)
      shared_ptr<Matrix> C(new Matrix((*U % *bvec + *bvec % *U)*0.5));

      // initial dR value.
      shared_ptr<Matrix> dR(new Matrix(*grad));
      shared_ptr<const Matrix> denom = compute_denom(C);
      *dR /= *denom;

      BFGS<Matrix> bfgs(denom);

      // solving Eq. 30 of Werner and Knowles
      // fixed U, solve for dR.
      for (int mmiter = 0; mmiter != max_mmicro_iter_; ++mmiter, ++tcount) {
        Timer mtimer;
        // first need to orthonomalize
        solver.orthog(dR);
        shared_ptr<Matrix> sigma = compute_sigma_R(jvec, dR, C, U);
        shared_ptr<Matrix> residual = solver.compute_residual(dR, sigma);

        const double error_mmicro = ::pow(residual->norm(),2.0) / residual->size();
        print_iteration(iter, miter, tcount, zero, error_mmicro, mtimer.tick());
        if (error_mmicro < thresh_mmicro_ || mmiter+1 == max_mmicro_iter_) { cout << endl; break; }

        // update dR;
        dR = bfgs.extrapolate(residual, solver.civec());
        dR->purify_redrotation(nclosed_,nact_,nvirt_);
      }

      // update U
      dR = solver.civec();
      shared_ptr<Matrix> dU = dR->exp();
      dU->purify_unitary();
      *U *= *dU;

    }

    U->purify_unitary();
    shared_ptr<const Coeff> newcc(new Coeff(*coeff_ * *U));
    coeff_ = newcc;

    resume_stdcout();
    print_iteration(iter, miter, tcount, energy_, error, timer.tick());
    mute_stdcout();

    if (error < thresh_) break;

  }
  resume_stdcout();

}




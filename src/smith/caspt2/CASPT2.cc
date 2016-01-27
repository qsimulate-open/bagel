//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <iostream>
#include <iomanip>
#include <src/smith/caspt2/CASPT2.h>
#include <src/util/math/linearRM.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

CASPT2::CASPT2::CASPT2(shared_ptr<const SMITH_Info<double>> ref) : SpinFreeMethod(ref) {
  eig_ = f1_->diag();
  nstates_ = ref->ciwfn()->nstates();
 
// from ss-caspt2 version for t2, r and s. 
//  t2 = init_amplitude(); 
//  r = init_residual();
//  s = init_residual();

// MS-CASPT2: t2 and s as MultiTensor (t2all, sall)
  for (int i = 0; i != nstates_; ++i) {
    auto tmp = make_shared<MultiTensor>(nstates_);
    for (auto& j : *tmp) 
      j = init_amplitude();
    t2all_.push_back(tmp);

   auto tmp2 = make_shared<MultiTensor>(nstates_);
    for (auto& j : *tmp2)
      j = init_residual();
    sall_.push_back(tmp2);
// added rall here.... maybe wrong
    rall_.push_back(tmp2->copy());
  }
}


void CASPT2::CASPT2::solve() {
  Timer timer;
  print_iteration();

// <proj| H |0> set to s in ss-caspt2 
  //shared_ptr<Queue> sourceq = make_sourceq();
  //while (!sourceq->done())
    //sourceq->next_compute();

//<proj_jst|H|0_K> set to sall in ms-caspt2 
  for (int istate = 0; istate != nstates_; ++istate) { //K states
    const double refen = info_->ciwfn()->energy(istate); 
    t2all_[istate]->fac(istate) = 0.0;
    sall_[istate]->fac(istate)  = 0.0;

    for (int jst=0; jst != nstates_; ++jst) { // <jst| 
      set_rdm(jst, istate);
      s = sall_[istate]->at(jst);
      shared_ptr<Queue> sourceq = make_sourceq(false, jst == istate);
      while(!sourceq->done())
        sourceq->next_compute();
    }
  } 

// set t2 to zero
 {
   for (int i = 0; i != nstates_; ++i) {
    t2all_[i]->zero();
    update_amplitude(t2all_[i], sall_[i]);
   }
  }

//maybe wrong... solver still not working
 LinearRM<Tensor> solver(30, s);

  Timer mtimer;
  int iter = 0;
  for ( ; iter != info_->maxiter(); ++iter) {
   
//ss-caspt2: R = <proj| H0 - E0 |1> + <proj | H |0> is set to r
//  first term
  //    shared_ptr<Queue> queue = make_residualq();
  //    while (!queue->done())
  //      queue->next_compute();
  //    diagonal(r, t2);

//ms-caspt2: R_K = <proj_jst| H0 - E0_K |1_ist> + <proj_jst| H |0_K> is set to rall
    //loop over state of interest
    for (int i = 0; i != nstates_; ++i) {  // K states

      //compute residuals named r for each K
      for (int ist = 0; ist != nstates_; ++ist) { // ist ket vector
        for (int jst = 0; jst != nstates_; ++jst) { // jst bra vector
        //first term <proj_jst| H0 - E0_K |1_ist>
           set_rdm(jst, ist);
           t2 = t2all_[i]->at(ist);
           r = rall_[i]->at(jst);
           shared_ptr<Queue> queue = make_residualq(false, jst == ist);
           while (!queue->done())
             queue->next_compute();
           diagonal(r, t2);

//getting errors with solver lines...
//           r = solver.compute_residual(t2, r);
//           t2 = solver.civec();

//cout << setprecision(10) << setw(20) << fixed << rall_[i]->rms() << "   rall in loop for state "<< i << endl;
//cout << setprecision(10) << setw(20) << fixed << t2all_[i]->rms() << "   t2all in loop for state  " << i << endl;

        }
      }
    }
  }


#if 0
    // E is now Hylleraas energy 
    // E = <1| H |0>
    // E += <1| H0 - E0 |1> + <1| H |0>
    for (int istate = 0; istate != nstates_; ++istate) {
      for (int ist = 0; ist != nstates_; ++ist) { // ket sector
        for (int jst = 0; jst != nstates_; ++jst) { // bra sector
          energy_ = detail::real(dot_product_transpose(sall_[istate], t2all_[istate]));
          cout << setprecision(10) << setw(20) << fixed << energy_[istate] << "E in loop " << endl;
          energy_[istate] += detail::real(dot_product_transpose(rall_[istate], t2all_[istate]));

    // get the root mean square
  for (int i = 0; i != nstates_; ++i) {
    const double err = rall_[i]->rms();
    print_iteration(iter, energy_[i], err, mtimer.tick());

    // computing delta T
    t2all_[i]->zero();
    update_amplitude(t2all_[i], rall_[i]);
    // zeroing out the residual
    rall_[i]->zero();
    if (err < info_->thresh()) break;
  }
  print_iteration(iter == info_->maxiter());
  timer.tick_print("CASPT2 energy evaluation");
//  cout << "    * CASPT2 energy : " << fixed << setw(20) << setprecision(10) << energy_+info_->ciwfn()->energy(0) << endl;
#endif
}

void CASPT2::CASPT2::solve_deriv() {

#if 0   
  Timer timer;
  shared_ptr<Queue> corrq = make_corrq();
  correlated_norm_ = accumulate(corrq);
  timer.tick_print("T1 norm evaluation");

  den2 = h1_->clone();
  shared_ptr<Queue> dens2 = make_densityq();
  while (!dens2->done())
    dens2->next_compute();

  den1 = h1_->clone();
  shared_ptr<Queue> dens1 = make_density1q();
  while (!dens1->done())
    dens1->next_compute();

  Den1 = init_residual();
  shared_ptr<Queue> Dens1 = make_density2q();
  while (!Dens1->done())
    Dens1->next_compute();
  timer.tick_print("Correlated density matrix evaluation");

  deci = make_shared<Tensor>(vector<IndexRange>{ci_});
  shared_ptr<Queue> dec = make_deciq();
  while (!dec->done())
    dec->next_compute();
  timer.tick_print("CI derivative evaluation");
  cout << endl;
#endif
 }


#endif

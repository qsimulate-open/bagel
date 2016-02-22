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
    rall_.push_back(tmp2->copy());
  }
}


void CASPT2::CASPT2::solve() {
  Timer timer;
  print_iteration();

  // <proj_jst|H|0_K> set to sall in ms-caspt2
  for (int istate = 0; istate != nstates_; ++istate) { //K states
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
      e0_ = e0all_[i];
      update_amplitude(t2all_[i], sall_[i]);
    }
  }

  // Linear solver for each state
  vector<shared_ptr<LinearRM<MultiTensor>>> solvers(nstates_);
  for (int i = 0; i != nstates_; ++i)
    solvers[i] = make_shared<LinearRM<MultiTensor>>(30, sall_[i]);

  energy_.resize(nstates_);
  vector<double> error(nstates_);
  Timer mtimer;
  int iter = 0;
  vector<bool> conv(nstates_, false);
  for ( ; iter != info_->maxiter(); ++iter) {
    // ms-caspt2: R_K = <proj_jst| H0 - E0_K |1_ist> + <proj_jst| H |0_K> is set to rall
    // loop over state of interest

    for (int i = 0; i != nstates_; ++i) {  // K states
      if (conv[i]) {
        print_iteration(iter, energy_[i], error[i], 0.0);
        solvers[i].reset();
        continue;
      }
      const double norm = t2all_[i]->norm();
      t2all_[i]->scale(1.0/norm);

      // compute residuals named r for each K
      for (int ist = 0; ist != nstates_; ++ist) { // ist ket vector
        for (int jst = 0; jst != nstates_; ++jst) { // jst bra vector
          // first term <proj_jst| H0 - E0_K |1_ist>
          set_rdm(jst, ist);
          t2 = t2all_[i]->at(ist);
          r = rall_[i]->at(jst);
          e0_ = e0all_[i];
          shared_ptr<Queue> queue = make_residualq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();
          diagonal(r, t2, jst == ist);
        }
      }

      // solve using subspace updates
      rall_[i] = solvers[i]->compute_residual(t2all_[i], rall_[i]);
      t2all_[i] = solvers[i]->civec();

      // energy is now the Hylleraas energy
      energy_[i] = detail::real(dot_product_transpose(sall_[i], t2all_[i]));
      energy_[i] += detail::real(dot_product_transpose(rall_[i], t2all_[i]));

      // compute rms for state i
      error[i] = rall_[i]->rms();
      print_iteration(iter, energy_[i], error[i], mtimer.tick());
      conv[i] = error[i] < info_->thresh();

      // compute delta t2 and update amplitude
      if (!conv[i]) {
        t2all_[i]->zero();
        update_amplitude(t2all_[i], rall_[i]);
        rall_[i]->zero();
      }
    }
    if (all_of(conv.begin(), conv.end(), [](bool i){ return i; })) break;
    if (nstates_ > 1) cout << endl;
  }
  print_iteration(iter == info_->maxiter());
  timer.tick_print("CASPT2 energy evaluation");
  cout << endl;
  for (int istate = 0; istate != nstates_; ++istate)
    cout << "    * CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << energy_[istate]+(*eref_)(istate,istate) << endl;

  // MS-CASPT2
  if (info_->do_ms() && nstates_ > 1) {
    Matrix fmn(nstates_, nstates_);

    for (int ist = 0; ist != nstates_; ++ist) {
      for (int jst = 0; jst != nstates_; ++jst) {
        if (ist == jst) {
          // set diagonal elements
          fmn(ist, ist) = energy_[ist] + (*eref_)(ist, ist);
        } else if (ist < jst) {
          // set off-diag elements
          // 1/2 [ <1g | H | Oe> + <0g |H | 1e > }
          fmn(jst, ist) = 0.5*(detail::real(dot_product_transpose(sall_[ist], t2all_[jst]))
                             + detail::real(dot_product_transpose(sall_[jst], t2all_[ist])))
                        + (*eref_)(jst, ist);
          fmn(ist, jst) = fmn(jst, ist);
        }
      }
    }

    // print out the effective Hamiltonian
    cout << endl;
    cout << "    * MS-CASPT2 Heff";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << fmn(ist, jst);
    }
    cout << endl << endl;

    VectorB eig(nstates_);
    fmn.diagonalize(eig);

    // energy printout
    for (int istate = 0; istate != nstates_; ++istate) {
      energy_[istate] = eig[istate];
      cout << "    * MS-CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << eig[istate] << endl;
    }
  }
}

void CASPT2::CASPT2::solve_deriv() {
  assert(nstates_ == 1);
  t2 = t2all_[0]->at(0);

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
  deci->allocate();
  shared_ptr<Queue> dec = make_deciq();
  while (!dec->done())
    dec->next_compute();
  timer.tick_print("CI derivative evaluation");
  cout << endl;
 }


#endif

////
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
    // TODO to be allocated in the gradient code
    lall_.push_back(tmp->copy());

    auto tmp2 = make_shared<MultiTensor>(nstates_);
    for (auto& j : *tmp2)
      j = init_residual();
    sall_.push_back(tmp2);
    rall_.push_back(tmp2->copy());
  }
  energy_.resize(nstates_);
  pt2energy_.resize(nstates_);
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

  // solve linear equation for t amplitudes
  t2all_ = solve_linear(sall_, t2all_);
  timer.tick_print("CASPT2 energy evaluation");
  cout << endl;

  for (int istate = 0; istate != nstates_; ++istate) {
    if (info_->shift() == 0) {
      pt2energy_[istate] = energy_[istate]+(*eref_)(istate,istate);
      cout << "    * CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] <<endl;
    } else {
      // will be used in normq
      n = init_residual();
      double norm = 0.0;
      for (int jst = 0; jst != nstates_; ++jst) { // bra
        for (int ist = 0; ist != nstates_; ++ist) { // ket
          set_rdm(jst, ist);
          t2 = t2all_[istate]->at(ist);
          shared_ptr<Queue> normq = make_normq(true, jst == ist);
          while (!normq->done())
            normq->next_compute();
          norm += dot_product_transpose(n, t2all_[istate]->at(jst));
        }
      }

      pt2energy_[istate] = energy_[istate]+(*eref_)(istate,istate) - info_->shift()*norm;
      cout << "    * CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
      cout << "        w/o shift correction  " << fixed << setw(20) << setprecision(10) << energy_[istate]+(*eref_)(istate,istate) <<endl;
      cout <<endl;
    }
  }

  // MS-CASPT2
// TODO move heff_ to the header ////////
shared_ptr<Matrix> heff_;
/////////////////////////////////////////
  if (info_->do_ms() && nstates_ > 1) {
    heff_ = make_shared<Matrix>(nstates_, nstates_);

    for (int ist = 0; ist != nstates_; ++ist) {
      for (int jst = 0; jst != nstates_; ++jst) {
        if (ist == jst) {
          // set diagonal elements
          (*heff_)(ist, ist) = pt2energy_[ist];
        } else if (ist < jst) {
          // set off-diag elements
          // 1/2 [ <1g | H | Oe> + <0g |H | 1e > }
          (*heff_)(jst, ist) = 0.5*(detail::real(dot_product_transpose(sall_[ist], t2all_[jst]))
                                  + detail::real(dot_product_transpose(sall_[jst], t2all_[ist])))
                             + (*eref_)(jst, ist);
          (*heff_)(ist, jst) = (*heff_)(jst, ist);
        }
      }
    }

    // print out the effective Hamiltonian
    cout << endl;
    cout << "    * MS-CASPT2 Heff";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    VectorB eig(nstates_);
    heff_->diagonalize(eig);

    // print out the eigen vector
    cout << endl;
    cout << "    * MS-CASPT2 rotation matrix";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    // energy printout
    for (int istate = 0; istate != nstates_; ++istate) {
      energy_[istate] = eig[istate];
      cout << "    * MS-CASPT2 energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << eig[istate] << endl;
    }
    cout << endl << endl;
  }

  // TODO
  // TODO move to the gradient code
  if (info_->do_ms() && nstates_ > 1) {
    // lamda eqn :   T_M <omega' | H | M' > T_M' + <omega' | f - E0_M + Eshift | Omega> lamdba  - E_shift * (T_M)^2 * <proj|Psi_M>= 0
    // compute first term and shift term (if used)
    print_iteration();

    // rall_[0] stores the result of summation over M'
// TODO
const int target = 0;
//  const int target = info_->target();
    rall_[0]->zero();
    for (int istate = 0; istate != nstates_; ++istate) //K states
      rall_[0]->ax_plus_y((*heff_)(istate, target), sall_[istate]);

    for (int istate = 0; istate != nstates_; ++istate) { //K states
      sall_[istate]->zero();
      sall_[istate]->ax_plus_y((*heff_)(istate, target), rall_[0]);
      if (info_->shift() != 0.0) {
        // subtract 2*Eshift*T_M^2*<proj|Psi_M> from source term
        n = init_residual();
        for (int jst = 0; jst != nstates_; ++jst) { // bra
          for (int ist = 0; ist != nstates_; ++ist) { // ket
            set_rdm(jst, ist);
            t2 = t2all_[istate]->at(ist);
            shared_ptr<Queue> normq = make_normq(true, jst == ist);
            while (!normq->done())
              normq->next_compute();

            sall_[istate]->at(jst)->ax_plus_y(-2.0 * info_->shift() * pow((*heff_)(istate, target), 2.0), n);
          }
        }
      }
    }

    // solve linear equation and store lambda in lall
    lall_ = solve_linear(sall_, lall_);
  }
}


// function to solve linear equation
vector<shared_ptr<MultiTensor_<double>>> CASPT2::CASPT2::solve_linear(vector<shared_ptr<MultiTensor_<double>>> s, vector<shared_ptr<MultiTensor_<double>>> t) {
  vector<shared_ptr<LinearRM<MultiTensor>>> solvers(nstates_);
  for (int i = 0; i != nstates_; ++i)
    solvers[i] = make_shared<LinearRM<MultiTensor>>(30, s[i]);

  // set t2 to guess vectors
  for (int i = 0; i != nstates_; ++i) {
    t[i]->zero();
    e0_ = e0all_[i] - info_->shift();
    update_amplitude(t[i], s[i]);
  }

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
      } else {
        rall_[i]->zero();
      }

      const double norm = t[i]->norm();
      t[i]->scale(1.0/norm);

      // compute residuals named r for each K
      for (int ist = 0; ist != nstates_; ++ist) { // ist ket vector
        for (int jst = 0; jst != nstates_; ++jst) { // jst bra vector
          // first term <proj_jst| H0 - E0_K |1_ist>
          set_rdm(jst, ist);
          t2 = t[i]->at(ist);
          r = rall_[i]->at(jst);

          e0_ = e0all_[i] - info_->shift();
          shared_ptr<Queue> queue = make_residualq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();
          diagonal(r, t2, jst == ist);
        }
      }
      // solve using subspace updates
      rall_[i] = solvers[i]->compute_residual(t[i], rall_[i]);
      t[i] = solvers[i]->civec();

      // energy is now the Hylleraas energy
      energy_[i] = detail::real(dot_product_transpose(s[i], t[i]));
      energy_[i] += detail::real(dot_product_transpose(rall_[i], t[i]));

      // compute rms for state i
      error[i] = rall_[i]->rms();
      print_iteration(iter, energy_[i], error[i], mtimer.tick());
      conv[i] = error[i] < info_->thresh();

      // compute delta t2 and update amplitude
      if (!conv[i]) {
        t[i]->zero();
        update_amplitude(t[i], rall_[i]);
      }
    }
    if (all_of(conv.begin(), conv.end(), [](bool i){ return i; })) break;
    if (nstates_ > 1) cout << endl;
  }
  print_iteration(iter == info_->maxiter());
  return t;
}


void CASPT2::CASPT2::solve_deriv() {
  assert(nstates_ == 1);
  t2 = t2all_[0]->at(0);

  Timer timer;
  {
    n = init_residual();
    shared_ptr<Queue> normq = make_normq();
    while (!normq->done())
      normq->next_compute();
    correlated_norm_ = dot_product_transpose(n, t2);
  }
  timer.tick_print("T1 norm evaluation");

  {
    den2 = h1_->clone();
    shared_ptr<Queue> dens2 = make_densityq();
    while (!dens2->done())
      dens2->next_compute();
  }

  {
    den1 = h1_->clone();
    shared_ptr<Queue> dens1 = make_density1q();
    while (!dens1->done())
      dens1->next_compute();
  }

  {
    Den1 = init_residual();
    shared_ptr<Queue> Dens1 = make_density2q();
    while (!Dens1->done())
      Dens1->next_compute();
  }
  timer.tick_print("Correlated density matrix evaluation");

  // rdm ci derivatives. Only for gradient computations
  // TODO this is only valid for single-state CASPT2
  {
    // first make ci_deriv_
    ci_deriv_ = make_shared<Civec>(info_->ref()->ciwfn()->det());
    const size_t cisize = ci_deriv_->size();

    const size_t cimax = 1000; // TODO interface to the input
    const size_t npass = (cisize-1)/cimax+1;
    const size_t chunk = (cisize-1)/npass+1;

    for (int ipass = 0; ipass != npass; ++ipass) {
      const size_t offset = ipass * chunk;
      const size_t size = min(chunk, cisize-offset);

      feed_rdm_deriv(offset, size); // this set ci_
      deci = make_shared<Tensor>(vector<IndexRange>{ci_});
      deci->allocate();
      shared_ptr<Queue> dec = make_deciq();
      while (!dec->done())
        dec->next_compute();
      copy_n(deci->vectorb()->data(), size, ci_deriv_->data()+offset);

      stringstream ss; ss << "CI derivative evaluation " << setw(5) << ipass+1 << " / " << npass;
      timer.tick_print(ss.str());
    }
    cout << endl;
  }
}

#endif

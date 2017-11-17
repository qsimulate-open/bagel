//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASA.cc
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
#include <src/smith/casa/CASA.h>
#include <src/util/math/linearRM.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

CASA::CASA::CASA(shared_ptr<const SMITH_Info<double>> ref) : SpinFreeMethod(ref) {
  eig_ = f1_->diag();
  nstates_ = ref->nact() ? ref->ciwfn()->nstates() : 1;

  // MS-CASA: t2 and s as MultiTensor (t2all, sall)
  for (int i = 0; i != nstates_; ++i) {
    auto tmp = make_shared<MultiTensor>(nstates_);
    auto tmp2 = make_shared<MultiTensor>(nstates_);
    for (int j = 0; j != nstates_; ++j)
      if (!info_->sssr() || i == j) {
        (*tmp)[j] = init_amplitude();
        (*tmp2)[j] = init_residual();
      }
    t2all_.push_back(tmp);
    sall_.push_back(tmp2);
    rall_.push_back(tmp2->clone());
  }
  energy_.resize(nstates_);
  pt2energy_.resize(nstates_);
}


void CASA::CASA::solve() {
  Timer timer;
  print_iteration();

  // <proj_jst|H|0_K> set to sall in ms-caspt2
  for (int istate = 0; istate != nstates_; ++istate) { //K states
    t2all_[istate]->fac(istate) = 0.0;
    sall_[istate]->fac(istate)  = 0.0;

    for (int jst=0; jst != nstates_; ++jst) { // <jst|
      if (info_->sssr() && jst != istate)
        continue;
      set_rdm(jst, istate);
      s = sall_[istate]->at(jst);
      shared_ptr<Queue> sourceq = make_sourceq(false, jst == istate);
      while(!sourceq->done())
        sourceq->next_compute();
    }
  }

  // solve linear equation for t amplitudes
  t2all_ = solve_linear(sall_, t2all_);
  timer.tick_print("CASA energy evaluation");
  cout << endl;

  // Compute and print energy
  for (int istate = 0; istate != nstates_; ++istate) {
    pt2energy_[istate] = energy_[istate]+(*eref_)(istate,istate);
    cout << "    * CASA energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] <<endl;
  }

  // MS-CASA
  if (info_->do_ms() && nstates_ > 1) {
    heff_ = make_shared<Matrix>(nstates_, nstates_);

    for (int ist = 0; ist != nstates_; ++ist) {
      auto sist = make_shared<MultiTensor>(nstates_);
      for (int jst = 0; jst != nstates_; ++jst) {
        if (sall_[ist]->at(jst)) {
          sist->at(jst) = sall_[ist]->at(jst);
        } else {
          set_rdm(jst, ist);
          s = init_residual();
          shared_ptr<Queue> sourceq = make_sourceq(false, jst == ist);
          while(!sourceq->done())
            sourceq->next_compute();
          sist->at(jst) = s;
        }
      }

      for (int jst = 0; jst != nstates_; ++jst) {
        if (ist == jst) {
          // set diagonal elements
          (*heff_)(ist, ist) = pt2energy_[ist];
        } else {
          // set off-diag elements
          // 1/2 [ <1g | H | Oe> + <0g |H | 1e > ]
          (*heff_)(jst, ist) = dot_product_transpose(sist, t2all_[jst]) + (*eref_)(jst, ist);
        }
      }
    }
    heff_->symmetrize();

    // print out the effective Hamiltonian
    cout << endl;
    cout << "    * MS-CASA Heff";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    VectorB eig(nstates_);
    heff_->diagonalize(eig);
    copy_n(eig.data(), nstates_, pt2energy_.data());

    // print out the eigen vector
    cout << endl;
    cout << "    * MS-CASA rotation matrix";
    for (int ist = 0; ist != nstates_; ++ist) {
      cout << endl << "      ";
      for (int jst = 0; jst != nstates_; ++jst)
        cout << setw(20) << setprecision(10) << (*heff_)(ist, jst);
    }
    cout << endl << endl;

    // energy printout
    for (int istate = 0; istate != nstates_; ++istate)
      cout << "    * MS-CASA energy : state " << setw(2) << istate << fixed << setw(20) << setprecision(10) << pt2energy_[istate] << endl;
    cout << endl << endl;
  } else {
    heff_ = make_shared<Matrix>(1,1);
    heff_->element(0,0) = 1.0;
  }
  energy_ = pt2energy_;
}


// function to solve linear equation
vector<shared_ptr<MultiTensor_<double>>> CASA::CASA::solve_linear(vector<shared_ptr<MultiTensor_<double>>> s, vector<shared_ptr<MultiTensor_<double>>> t) {
  Timer mtimer;
  // ms-caspt2: R_K = <proj_jst| H0 - E0_K |1_ist> + <proj_jst| H |0_K> is set to rall
  // loop over state of interest
  bool converged = true;

  e0f_ = e0all_;
  for (int ist = 0; ist != nstates_; ++ist)
    e0all_[ist] = info_->ref()->energy(ist);

  for (int i = 0; i != nstates_; ++i) {  // K states
    bool conv = false;
    double error = 0.0;
    e0_ = e0all_[i] - info_->shift();
    energy_[i] = 0.0;
    // set guess vector
    t[i]->zero();
    if (s[i]->rms() < 1.0e-15) {
      print_iteration(0, 0.0, 0.0, mtimer.tick());
      if (i+1 != nstates_) cout << endl;
      continue;
    } else {
      update_amplitude_casa(t[i], s[i], i);
    }

    auto solver = make_shared<LinearRM<MultiTensor>>(info_->davidson_subspace(), s[i]);
    int iter = 0;
    for ( ; iter != info_->maxiter(); ++iter) {
      rall_[i]->zero();

      const double norm = t[i]->norm();
      t[i]->scale(1.0/norm);

      // compute residuals named r for each K
      for (int jst = 0; jst != nstates_; ++jst) { // jst bra vector
        for (int ist = 0; ist != nstates_; ++ist) { // ist ket vector
          if (info_->sssr() && (jst != i || ist != i))
            continue;
          // first term <proj_jst| H0 - E0_K |1_ist>
          set_rdm(jst, ist);
          t2 = t[i]->at(ist);
          r = rall_[i]->at(jst);

          // compute residuals named r for each K
          const double e0save = e0_;
          e0_ = e0all_[i] - e0all_[ist];
          shared_ptr<Queue> queue = make_residualq(false, jst == ist);
          while (!queue->done())
            queue->next_compute();
          diagonal(r, t2, jst == ist);
          e0_ = e0save;
        }
      }
      // solve using subspace updates
      rall_[i] = solver->compute_residual(t[i], rall_[i]);
      t[i] = solver->civec();

      // energy is now the Hylleraas energy
      energy_[i] = detail::real(dot_product_transpose(s[i], t[i]));
      energy_[i] += detail::real(dot_product_transpose(rall_[i], t[i]));

      // compute rms for state i
      error = rall_[i]->norm() / pow(rall_[i]->size(), 0.25);
      print_iteration(iter, energy_[i], error, mtimer.tick());
      conv = error < info_->thresh();

      // compute delta t2 and update amplitude
      if (!conv) {
        t[i]->zero();
        update_amplitude_casa(t[i], rall_[i], i);
      }
      if (conv) break;
    }
    if (i+1 != nstates_) cout << endl;
    converged &= conv;
  }
  print_iteration(!converged);
  return t;
}


// Update_amplitude is hard-coded to use the Fock matrix, so we substitute in CASPT2 e0 here
void CASA::CASA::update_amplitude_casa(std::shared_ptr<MultiTensor> t, std::shared_ptr<const MultiTensor> r, const int istate) {
  vector<double> e0allsave = e0all_;
  double e0save = e0_;

  e0all_ = e0f_;
  e0_ = e0all_[istate] - info_->shift();

  update_amplitude(t, r);

  e0all_ = e0allsave;
  e0_ = e0save;
}

void CASA::CASA::solve_gradient(const int targetJ, const int targetI, shared_ptr<const NacmType> nacmtype, const bool nocider) {
  throw runtime_error("Gradients are not implemented for CAS/A.");
}

#endif

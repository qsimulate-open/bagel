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


#include <src/smith/caspt2/CASPT2.h>
#include <src/util/math/linearRM.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

CASPT2::CASPT2::CASPT2(shared_ptr<const SMITH_Info<double>> ref) : SpinFreeMethod(ref) {
  eig_ = f1_->diag();
  t2 = init_amplitude();
  r = init_residual();
  s = init_residual();
}


void CASPT2::CASPT2::solve() {
  Timer timer;
  print_iteration();
  // <proj| H |0> set to s
  shared_ptr<Queue> sourceq = make_sourceq();
  while (!sourceq->done())
    sourceq->next_compute();
  LinearRM<Tensor> solver(30, s);

  Timer mtimer;
  int iter = 0;
  for ( ; iter != info_->maxiter(); ++iter) {
    // E = <1| H |0>
    energy_ = detail::real(dot_product_transpose(s, t2));

    // R = <proj| H0 - E0 |1> + <proj | H |0> is set to r
    //  first term
    shared_ptr<Queue> queue = make_residualq();
    while (!queue->done())
      queue->next_compute();
    diagonal(r, t2);

    // TODO this will be replaced by subspace updates
    //  update_amplitude(t2, r);
    solver.compute_residual(t2->copy(), s->copy());

    //  second term
    r->ax_plus_y(1.0, s);
    // E is now Hylleraas energy 
    // E += <1| H0 - E0 |1> + <1| H |0>
    energy_ += detail::real(dot_product_transpose(r, t2));

    // get the root mean square
    const double err = r->rms();
    print_iteration(iter, energy_, err, mtimer.tick());

    // zeroing out the residual
    r->zero();
    if (err < info_->thresh()) break;
  }
  print_iteration(iter == info_->maxiter());
  timer.tick_print("CASPT2 energy evaluation");
  cout << "    * CASPT2 energy : " << fixed << setw(20) << setprecision(10) << energy_+info_->ciwfn()->energy(0) << endl;
}


void CASPT2::CASPT2::solve_deriv() {
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
}


#endif

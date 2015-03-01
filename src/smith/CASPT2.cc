//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/CASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

CASPT2::CASPT2::CASPT2(shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {
  this->eig_ = f1_->diag();
  t2 = init_amplitude();
  e0_ = this->e0();
  r = t2->clone();
  den1 = h1_->clone();
  den2 = h1_->clone();
  Den1 = v2_->clone();
  deci = make_shared<Tensor>(vector<IndexRange>{ci_});
}


void CASPT2::CASPT2::solve() {
  Timer timer;
  this->print_iteration();
  int iter = 0;
  for ( ; iter != ref_->maxiter(); ++iter) {
    shared_ptr<Queue> energyq = make_energyq();
    this->energy_ = accumulate(energyq);
    shared_ptr<Queue> queue = make_residualq();
    while (!queue->done())
      queue->next_compute();
    diagonal(r, t2);
    this->energy_ += dot_product_transpose(r, t2) * 0.25;
    const double err = r->rms();
    this->print_iteration(iter, this->energy_, err);

    this->update_amplitude(t2, r);
    r->zero();
    if (err < ref_->thresh()) break;
  }
  this->print_iteration(iter == ref_->maxiter());
  timer.tick_print("CASPT2 energy evaluation");
}


void CASPT2::CASPT2::solve_deriv() {
  Timer timer;
  shared_ptr<Queue> corrq = make_corrq();
  correlated_norm_ = accumulate(corrq);
  timer.tick_print("T1 norm evaluation");

  shared_ptr<Queue> dens2 = make_densityq();
  while (!dens2->done())
    dens2->next_compute();
  shared_ptr<Queue> dens1 = make_density1q();
  while (!dens1->done())
    dens1->next_compute();
  shared_ptr<Queue> Dens1 = make_density2q();
  while (!Dens1->done())
    Dens1->next_compute();
  timer.tick_print("Correlated density matrix evaluation");

  shared_ptr<Queue> dec = make_deciq();
  while (!dec->done())
    dec->next_compute();
  timer.tick_print("CI derivative evaluation");
  cout << endl;
}


#endif

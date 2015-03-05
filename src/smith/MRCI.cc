//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI.cc
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


#include <src/util/math/davidson.h>
#include <src/smith/extrap.h>
#include <src/smith/MRCI.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

MRCI::MRCI::MRCI(shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {
  this->eig_ = f1_->diag();
  t2 = init_amplitude();
  r = t2->clone();
  s = t2->clone();
  n = t2->clone();
  den1 = h1_->clone();
  den2 = h1_->clone();
  Den1 = v2_->clone();
  deci = make_shared<Tensor>(vector<IndexRange>{ci_});
}


void MRCI::MRCI::solve() {
  Timer timer;
  this->print_iteration();

  auto queue = make_sourceq();
  while (!queue->done())
    queue->next_compute();
  queue = make_normq();
  while (!queue->done())
    queue->next_compute();

  DavidsonDiag_<Amplitude, Residual> davidson(1, 100);

  const double refen = ref_->ciwfn()->energy(ref_->target()) - this->core_energy_ - ref_->geom()->nuclear_repulsion();
  auto a0 = make_shared<Amplitude>(1.0, t2, n, this);
  auto r0 = make_shared<Residual>(refen, s, this);
  davidson.compute(a0, r0);
  r = davidson.residual().front()->tensor();
  this->update_amplitude(t2, r);

  int iter = 0;
  for ( ; iter != ref_->maxiter(); ++iter) {
    queue = make_residualq();
    while (!queue->done())
      queue->next_compute();

    queue = make_normq();
    while (!queue->done())
      queue->next_compute();

    a0 = make_shared<Amplitude>(0.0, t2, n, this);
    r0 = make_shared<Residual>(dot_product_transpose(s, t2), r, this);

    this->energy_ = davidson.compute(a0, r0);
    r = davidson.residual()[0]->tensor();
    const double err = r->rms();
    this->print_iteration(iter, this->energy_, err);

    t2->zero();
    this->update_amplitude(t2, r);

    if (err < ref_->thresh()) break;
  }
  this->print_iteration(iter == ref_->maxiter());
  timer.tick_print("MRCI energy evaluation");
}


void MRCI::MRCI::solve_deriv() {
  throw std::logic_error("Nuclear gradients not implemented for MRCI");
}


#endif

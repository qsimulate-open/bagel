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


#include <src/smith/MRCI.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

MRCI::MRCI::MRCI(shared_ptr<const SMITH_Info> ref) : SpinFreeMethod(ref) {
  this->eig_ = f1_->diag();
  t2 = init_amplitude();
  r = t2->clone();
  s = t2->clone();
  den1 = h1_->clone();
  den2 = h1_->clone();
  Den1 = v2_->clone();
  deci = make_shared<Tensor>(vector<IndexRange>{ci_});
}


void MRCI::MRCI::solve() {
  Timer timer;
  this->print_iteration();
  shared_ptr<Queue> sq = make_sourceq();
  while (!sq->done())
    sq->next_compute();
  const double ee = e0_;
  const double refen = ref_->ciwfn()->energy(ref_->target()) - this->core_energy_ - ref_->geom()->nuclear_repulsion();

  int iter = 0;
  for ( ; iter != ref_->maxiter(); ++iter) {
    e0_ = 0.0;
    shared_ptr<Queue> queue = make_residualq();
    while (!queue->done())
      queue->next_compute();
    r->ax_plus_y(2.0, s);
    this->energy_ = dot_product_transpose(r, t2) + refen;

    // norm
    shared_ptr<Queue> corrq = make_corrq();
    this->energy_ /= (1.0+accumulate(corrq));

    // compute residual
    e0_ = this->energy_;
    queue = make_residualq();
    while (!queue->done())
      queue->next_compute();
    r->ax_plus_y(1.0, s);

    const double err = r->rms();
    this->energy_ += this->core_energy_ + ref_->geom()->nuclear_repulsion();
    this->print_iteration(iter, this->energy_, err);

    e0_ = ee;
    this->update_amplitude(t2, r);
    r->zero();
    if (err < ref_->thresh()) break;
  }
  this->print_iteration(iter == ref_->maxiter());
  timer.tick_print("MRCI energy evaluation");
}


void MRCI::MRCI::solve_deriv() {
  throw std::logic_error("Nuclear gradients not implemented for MRCI");
}


#endif

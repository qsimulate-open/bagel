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
  e0_ = this->e0();
  r = t2->clone();
  den1 = h1_->clone();
  den2 = h1_->clone();
  Den1 = v2_->clone();
  deci = make_shared<Tensor>(vector<IndexRange>{ci_});
}


void MRCI::MRCI::solve() {
  Timer timer;
  this->print_iteration();
  int iter = 0;
  for ( ; iter != ref_->maxiter(); ++iter) {
    shared_ptr<Queue> queue = make_residualq();
    while (!queue->done())
      queue->next_compute();
    const double err = r->rms();
    this->print_iteration(iter, this->energy_, err);

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

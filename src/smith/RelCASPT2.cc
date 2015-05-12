//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2.cc
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


#include <src/smith/RelCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

RelCASPT2::RelCASPT2::RelCASPT2(shared_ptr<const SMITH_Info<std::complex<double>>> ref) : SpinFreeMethod(ref) {
  auto eig = f1_->diag();
  eig_.resize(eig.size());
  for (int i = 0; i != eig.size(); ++i)
    eig_[i] = real(eig[i]);
  t2 = init_amplitude();
  r = t2->clone();
}


void RelCASPT2::RelCASPT2::solve() {
  Timer timer;
  print_iteration();
  Timer mtimer;
  int iter = 0;
  for ( ; iter != info_->maxiter(); ++iter) {
    shared_ptr<Queue> energyq = make_energyq();
    energy_ = accumulate(energyq);
    shared_ptr<Queue> queue = make_residualq();
    while (!queue->done())
      queue->next_compute();
    energy_ += detail::real(dot_product_transpose(r, t2));
    const double err = r->rms();
    print_iteration(iter, energy_, err, mtimer.tick());

    update_amplitude(t2, r);
    r->zero();
    if (err < info_->thresh()) break;
  }
  print_iteration(iter == info_->maxiter());
  timer.tick_print("CASPT2 energy evaluation");
}


void RelCASPT2::RelCASPT2::solve_deriv() {
  throw std::logic_error("Nuclear gradients not implemented for RelCASPT2");
}


#endif

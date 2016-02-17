//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASPT2.cc
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


#include <src/smith/relcaspt2/RelCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

RelCASPT2::RelCASPT2::RelCASPT2(shared_ptr<const SMITH_Info<std::complex<double>>> ref) : SpinFreeMethod(ref) {
  auto eig = f1_->diag();
  eig_.resize(eig.size());
  for (int i = 0; i != eig.size(); ++i)
    eig_[i] = real(eig[i]);
  t2 = init_amplitude();
  r = init_residual();
  s = init_residual();
}


void RelCASPT2::RelCASPT2::solve() {
  Timer timer;
  print_iteration();
  shared_ptr<Queue> sourceq = make_sourceq();
  while (!sourceq->done())
    sourceq->next_compute();
  Timer mtimer;
  int iter = 0;
  energy_.resize(1);
  for ( ; iter != info_->maxiter(); ++iter) {
    energy_[0] = detail::real(dot_product_transpose(s, t2));
    shared_ptr<Queue> queue = make_residualq();
    while (!queue->done())
      queue->next_compute();
    diagonal(r, t2);
    r->ax_plus_y(1.0, s);
    energy_[0] += detail::real(dot_product_transpose(r, t2));
    const double err = r->rms();
    print_iteration(iter, energy_[0], err, mtimer.tick());

    update_amplitude(t2, r);
    r->zero();
    if (err < info_->thresh()) break;
  }
  print_iteration(iter == info_->maxiter());
  timer.tick_print("CASPT2 energy evaluation");
  cout << "    * CASPT2 energy : " << fixed << setw(20) << setprecision(10) << energy_[0]+info_->ciwfn()->energy(0) << endl;
}


void RelCASPT2::RelCASPT2::solve_deriv() {
  throw std::logic_error("Nuclear gradients not implemented for RelCASPT2");
}


#endif

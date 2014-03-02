//
// BAGEL - Parallel electron correlation program.
// Filename: harrison_denom.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <iomanip>
#include <stdexcept>
#include <bitset>
#include <src/fci/harrison.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>
#include <src/util/taskqueue.h>
#include <iostream>
#include <sstream>
#include <src/fci/hzdenomtask.h>

using namespace std;
using namespace bagel;


void HarrisonZarrabian::const_denom() {
  Timer denom_t;
  auto h = make_shared<Matrix>(norb_, 1);
  auto jop = make_shared<Matrix>(norb_, norb_);
  auto kop = make_shared<Matrix>(norb_, norb_);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop->element(j, i) = jop->element(i, j) = 0.5*jop_->mo2e_hz(j, i, j, i);
      kop->element(j, i) = kop->element(i, j) = 0.5*jop_->mo2e_hz(j, i, i, j);
    }
    h->element(i,0) = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<Civec>(det());

  double* iter = denom_->data();
  TaskQueue<HZDenomTask> tasks(det()->string_bits_a().size());
  for (auto& ia : det()->string_bits_a()) {
    tasks.emplace_back(iter, ia, det_, jop, kop, h);
    iter += det()->string_bits_b().size();
  }

  tasks.compute();
  denom_t.tick_print("denom");
}

void HarrisonZarrabian::update(shared_ptr<const Coeff> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, c, "HZ");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}

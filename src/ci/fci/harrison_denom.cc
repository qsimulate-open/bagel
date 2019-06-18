//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: harrison_denom.cc
// Copyright (C) 2011 Toru Shiozaki
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

#include <iomanip>
#include <stdexcept>
#include <bitset>
#include <src/ci/fci/harrison.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>
#include <src/util/taskqueue.h>
#include <iostream>
#include <sstream>
#include <src/ci/fci/hzdenomtask.h>

using namespace std;
using namespace bagel;


void HarrisonZarrabian::const_denom() {
  Timer denom_t;
  auto h = make_shared<VectorB>(norb_);
  auto jop = make_shared<Matrix>(norb_, norb_);
  auto kop = make_shared<Matrix>(norb_, norb_);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      (*jop)(j, i) = (*jop)(i, j) = 0.5*jop_->mo2e_hz(j, i, j, i);
      (*kop)(j, i) = (*kop)(i, j) = 0.5*jop_->mo2e_hz(j, i, i, j);
    }
    (*h)(i) = jop_->mo1e(i,i);
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

void HarrisonZarrabian::update(shared_ptr<const Matrix> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  coeff_ = c; 
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, coeff_, store_half_ints_, "HZ");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}

//
// BAGEL - Parallel electron correlation program.
// Filename: fci_init.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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
#include <src/fci/mofile.h>
#include <src/rysint/eribatch.h>
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
  unique_ptr<double[]> h(new double[norb_]);
  unique_ptr<double[]> jop(new double[norb_*norb_]);
  unique_ptr<double[]> kop(new double[norb_*norb_]);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, j, i);
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, i, j);
    }
    h[i] = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<Civec>(det());

  double* iter = denom_->data();
  vector<HZDenomTask> tasks;
  tasks.reserve(det()->stringa().size());
  for (auto& ia : det()->stringa()) {
    tasks.push_back(HZDenomTask(iter, ia, det_, jop.get(), kop.get(), h.get()));
    iter += det()->stringb().size();
  }

  TaskQueue<HZDenomTask> tq(tasks);
  tq.compute(resources__->max_num_threads());
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

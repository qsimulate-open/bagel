//
// BAGEL - Parallel electron correlation program.
// Filename: zharrison_denom.cc
// Copyright (C) 2013 Toru Shiozaki
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
#include <src/zfci/zharrison.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace bagel;


void ZHarrison::const_denom() {
  Timer denom_t;
  auto h = make_shared<ZMatrix>(norb_, 1);
  auto jaa = make_shared<ZMatrix>(norb_, norb_);
  auto jab = make_shared<ZMatrix>(norb_, norb_);
  auto jbb = make_shared<ZMatrix>(norb_, norb_);
  auto kaa = make_shared<ZMatrix>(norb_, norb_);
  auto kbb = make_shared<ZMatrix>(norb_, norb_);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j != norb_; ++j) {
      jaa->element(j, i) = 0.5*jop_->mo2e(bitset<4>("0000"), j, i, j, i);
      jab->element(j, i) = 0.5*jop_->mo2e(bitset<4>("0101"), j, i, j, i);
      jbb->element(j, i) = 0.5*jop_->mo2e(bitset<4>("1111"), j, i, j, i);
      kaa->element(j, i) = 0.5*jop_->mo2e(bitset<4>("0000"), j, i, i, j);
      kbb->element(j, i) = 0.5*jop_->mo2e(bitset<4>("1111"), j, i, i, j);
    }
    h->data(i) = jop_->mo1e(bitset<2>("00"), i,i);
    assert(abs(h->data(i) - jop_->mo1e(bitset<2>("11"), i,i)) < 1.0e-8);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<RelDvec>(space_, 1);

  const size_t est = accumulate(space_->detmap().begin(), space_->detmap().end(), 0ull, [](size_t r, pair<int,shared_ptr<Determinants>> i){ return r+i.second->stringa().size(); });
//TaskQueue<ZHarrisonDenomTask> tasks(est);

  for (auto& i : space_->detmap()) {
    shared_ptr<const Determinants> det = i.second; 
#if 0
    double* iter = denom_->data();
    for (auto& ia : det()->stringa()) {
      tasks.emplace_back(iter, ia, det_, jop.get(), kop.get(), h.get());
      iter += det()->stringb().size();
    }
#endif
  }

//tasks.compute();
  denom_t.tick_print("denom");
}

#if 0
void HarrisonZarrabian::update(shared_ptr<const Coeff> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, c, "HZ");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}
#endif

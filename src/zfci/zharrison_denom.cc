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
#include <src/zfci/zharrison.h>
#include <src/fci/hzdenomtask.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>
#include <iostream>
#include <sstream>

using namespace std;
using namespace bagel;


void ZHarrison::const_denom() {
  Timer denom_t;
  auto h = make_shared<Matrix>(norb_, 1);
  auto jop = make_shared<Matrix>(norb_, norb_);
  auto kop = make_shared<Matrix>(norb_, norb_);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j != norb_; ++j) {
      jop->element(j, i) = 0.5*jop_->mo2e("0000", j, i, j, i).real();
      kop->element(j, i) = 0.5*jop_->mo2e("1111", j, i, i, j).real();
      // assert for Kramers and symmetry
      assert(fabs(jop_->mo2e("0000", j, i, j, i).imag()) < 1.0e-8);
      assert(fabs(jop_->mo2e("1111", j, i, i, j).imag()) < 1.0e-8);
      assert(fabs(jop_->mo2e("0101", j, i, j, i).imag()) < 1.0e-8);
    }
    h->element(i, 0) = jop_->mo1e("00", i,i).real();
    // assert for Kramers and symmetry
    assert(abs(jop_->mo1e("00", i,i) - jop_->mo1e("11", i,i)) < 1.0e-8);
    assert(abs(jop_->mo1e("00", i,i).imag()) < 1.0e-8);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<RelDvec>(space_, 1);

  const size_t est = accumulate(space_->detmap().begin(), space_->detmap().end(), 0ull, [](size_t r, pair<int,shared_ptr<Determinants>> i){ return r+i.second->stringa().size(); });
  TaskQueue<HZDenomTask> tasks(est);

  for (auto& i : space_->detmap()) {
    shared_ptr<const Determinants> det = i.second;
    shared_ptr<Dvec> cdenom = denom_->find(det->nelea(), det->neleb());
    double* dptr = cdenom->data();
    for (auto& ia : det->stringa()) {
      tasks.emplace_back(dptr, ia, det, jop, kop, h);
      dptr += det->stringb().size();
    }
  }

  tasks.compute();
  denom_t.tick_print("denom");
}

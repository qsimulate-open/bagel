//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zharrison_denom.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/ci/zfci/zharrison.h>
#include <src/ci/fci/hzdenomtask.h>

using namespace std;
using namespace bagel;

void ZHarrison::const_denom() {
  Timer denom_t;
  auto h = make_shared<VectorB>(norb_);
  auto jop = make_shared<Matrix>(norb_, norb_);
  auto kop = make_shared<Matrix>(norb_, norb_);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j != norb_; ++j) {
      jop->element(j, i) = 0.25*jop_->mo2e("0000", j, i, j, i).real()
                         + 0.25*jop_->mo2e("1111", j, i, j, i).real();
      kop->element(j, i) = 0.25*jop_->mo2e("0000", j, i, i, j).real()
                         + 0.25*jop_->mo2e("1111", j, i, i, j).real();
    }
    (*h)(i) = 0.5*jop_->mo1e("00", i,i).real()
            + 0.5*jop_->mo1e("11", i,i).real();
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<RelDvec>(space_, 1);

  const size_t est = accumulate(space_->detmap().begin(), space_->detmap().end(), 0ull,
                                [](size_t r, pair<pair<int,int>,shared_ptr<Determinants>> i){ return r+i.second->string_bits_a().size(); });
  TaskQueue<HZDenomTask> tasks(est);

  for (auto& i : space_->detmap()) {
    shared_ptr<const Determinants> det = i.second;
    shared_ptr<Dvec> cdenom = denom_->find(det->nelea(), det->neleb());
    double* dptr = cdenom->data();
    for (auto& ia : det->string_bits_a()) {
      tasks.emplace_back(dptr, ia, det, jop, kop, h);
      dptr += det->string_bits_b().size();
    }
  }

  tasks.compute();
  denom_t.tick_print("denom");
}

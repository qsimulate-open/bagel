//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gradbatch.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <src/integral/rys/gradbatch.h>
#include <src/integral/rys/inline.h>
#include <src/integral/rys/erirootlist.h>
#include <src/integral/rys/breitrootlist.h>
#include <src/integral/rys/spin2rootlist.h>


using namespace std;
using namespace bagel;

GradBatch::GradBatch(const array<shared_ptr<const Shell>,4>& shells, const double max_density, const double dummy, const bool dum,
                     shared_ptr<StackMem> stack)
 : ERIBatch_base(shells, 1, 0, stack) {

  const double integral_thresh = (max_density != 0.0) ? (PRIM_SCREEN_THRESH / max_density) : 0.0;
  compute_ssss(integral_thresh);

  centers_ = 4;
  for (auto& i : shells) if (i->dummy()) --centers_;

  set_exponents();

  root_weight(this->primsize_);
}


void GradBatch::set_exponents() {
  exponents_ = unique_ptr<double[]>(new double[primsize_*4]);
  double* tmp = exponents_.get();
  for (auto& i0 : basisinfo_[0]->exponents()) {
    for (auto& i1 : basisinfo_[1]->exponents()) {
      for (auto& i2 : basisinfo_[2]->exponents()) {
        for (auto& i3 : basisinfo_[3]->exponents()) {
          tmp[0] = i0;
          tmp[1] = i1;
          tmp[2] = i2;
          tmp[3] = i3;
          tmp += 4;
        }
      }
    }
  }
}

void GradBatch::root_weight(const int ps) {
  if (amax_ + cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int i = screening_[j];
      if (T_[i] < 1.0e-8) {
        weights_[i] = 1.0;
      } else {
        const double sqrtt = sqrt(T_[i]);
        const double erfsqt = inline_erf(sqrtt);
        weights_[i] = erfsqt * sqrt(pi__) * 0.5 / sqrtt;
      }
    }
  } else {
    eriroot__.root(rank_, T_, roots_, weights_, ps);
  }
}

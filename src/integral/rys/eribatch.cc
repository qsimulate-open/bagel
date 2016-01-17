//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: eribatch.cc
// Copyright (C) 2009 Toru Shiozaki
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

#include <src/integral/rys/eribatch.h>
#include <src/integral/rys/inline.h>
#include <src/integral/rys/erirootlist.h>

using namespace std;
using namespace bagel;

ERIBatch::ERIBatch(const array<shared_ptr<const Shell>,4>& _info, const double max_density, const double dummy, const bool dum,
                   shared_ptr<StackMem> stack) :  ERIBatch_base(_info, 0, 0, stack) {

  const double integral_thresh = (max_density != 0.0) ? (PRIM_SCREEN_THRESH / max_density) : 0.0;
  compute_ssss(integral_thresh);

#ifdef LIBINT_INTERFACE
  assert(false);
#endif

  root_weight(this->primsize_);
}

void ERIBatch::root_weight(const int ps) {
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


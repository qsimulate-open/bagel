//
// BAGEL - Parallel electron correlation program.
// Filename: gradbatch.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <src/integral/rys/gradbatch.h>
#include <src/integral/rys/inline.h>

using namespace std;
using namespace bagel;

GradBatch::GradBatch(const array<shared_ptr<const Shell>,4>& shells, const double max_density, const double dummy, const bool dum,
                     shared_ptr<StackMem> stack)
 : ERIBatch_base(shells, max_density, 1, 0, stack) {

  centers_ = 4;
  for (auto& i : shells) if (i->dummy()) --centers_;

  set_exponents();
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

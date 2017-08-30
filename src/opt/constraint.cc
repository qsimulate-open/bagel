//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: constraint.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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

#include <src/opt/constraint.h>

using namespace std;
using namespace bagel;

OptExpBonds::OptExpBonds(shared_ptr<const PTree> inp) {
  pair_ = inp->get_array<int,2>("pair");
  if ((pair_[0] > pair_[1])) {
    int tmp = pair_[1];
    pair_[1] = pair_[0];
    pair_[0] = tmp;
  }
  for (int p = 0; p != 2; ++p) pair_[p]--;

  cout << "  * Explicit bond pair added between " << pair_[0] << " and " << pair_[1] << endl;
}

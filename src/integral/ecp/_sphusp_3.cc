//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _sphusp_3.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#include <algorithm>
#include <cassert>
#include <src/integral/ecp/sphusplist.h>

using namespace std;
using namespace bagel;

vector<double> SphUSPList::sphusp_3(const int m) {

  vector<double> c;
  constexpr double coeff[70] = {   0.000000000000000e+00,   1.770130769779930e+00,   0.000000000000000e+00,  -5.900435899266435e-01,   0.000000000000000e+00,
   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,
   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,
   2.890611442640554e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,
   0.000000000000000e+00,  -4.570457994644658e-01,   0.000000000000000e+00,  -4.570457994644658e-01,   0.000000000000000e+00,
   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   1.828183197857863e+00,   0.000000000000000e+00,
   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,  -1.119528997770346e+00,
   0.000000000000000e+00,  -1.119528997770346e+00,   0.000000000000000e+00,   0.000000000000000e+00,   7.463526651802308e-01,
  -4.570457994644658e-01,   0.000000000000000e+00,  -4.570457994644658e-01,   0.000000000000000e+00,   0.000000000000000e+00,
   0.000000000000000e+00,   0.000000000000000e+00,   1.828183197857863e+00,   0.000000000000000e+00,   0.000000000000000e+00,
   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   1.445305721320277e+00,
   0.000000000000000e+00,  -1.445305721320277e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,
   5.900435899266435e-01,   0.000000000000000e+00,  -1.770130769779930e+00,   0.000000000000000e+00,   0.000000000000000e+00,
   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00,   0.000000000000000e+00
};

  assert(abs(m) <= 3);
  const int size_c = (3 + 1) * (3 + 2) / 2;
  const int mu = m + 3;
  const int i0 = mu * size_c;
  for (int i = i0; i != i0 + size_c; ++i) c.push_back(coeff[i]);
  return c;

}

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: molden_transforms.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

//
// This was originally code found in src/integral/.../carsph_gen.cc
// but repurposed here for the molden section. At the moment, it is limited to at most
// f functions.
//

#include <src/util/math/comb.h>
#include <src/util/math/factorial.h>
#include <src/util/io/moldenin.h>

using namespace std;
using namespace bagel;

#define LARGE 32
#define LEND 5

const static Comb comb;
const static Factorial factorial;

void MoldenIn::compute_transforms() {
  vector<pair<int, double>> s0(1, make_pair(0, 1.0));
  vector<vector<pair<int, double>>> s1(1, s0);
  lmtuv_.push_back(s1);

  for (int l = 1; l != LEND; ++l) {

    map<int, int> mapping;
    int cnt = 0;
    for (int z = 0; z <= l; ++z) {
      for (int y = 0; y <= l - z; ++y) {
        const int x = l - y - z;
        if (x < 0) continue;
        const int key = x + y * LARGE + z * LARGE * LARGE;
        mapping.emplace(key, cnt);
        ++cnt;
      }
    }

    vector<vector<pair<int, double>>> mtuv;
    for (int n = 0; n != 2 * l + 1; ++n) {
      int m = l - (n / 2);
      if (n % 2 == 1) m *= -1;

      const int vm2 = m < 0 ? 1 : 0;
      const int absm = m > 0 ? m : -m;
      vector<pair<int, double>> tuv;

      const double Nlms = 1.0 / pow(2.0, absm) / factorial(l) * sqrt((m == 0 ? 1.0 : 2.0) * factorial(l + absm) * factorial(l - absm));
      const int tmax = (l - absm) / 2;
      for (int t = 0; t <= tmax; ++t) {
        for (int u = 0; u <= t; ++u) {
          const int vmax = 2 * ((absm - vm2) / 2) + vm2;
          for (int v2 = vm2; v2 <= vmax; v2 += 2) {
            assert((v2 - vm2) % 2 == 0);
            const double Clmtuv = pow(-1.0, t + ((v2 - vm2) / 2)) * pow(0.25, t)
                                * comb(l, t) * comb(l - t, absm + t)
                                * comb(t, u) * comb(absm, v2) * Nlms;
            const int xexp = 2 * t + absm - 2 * u - v2;
            const int yexp = 2 * u + v2;
            const int zexp = l - 2 * t - absm;
            auto current = mapping.find(xexp + yexp * LARGE + zexp * LARGE * LARGE);
            assert(current != mapping.end());
            tuv.push_back({current->second, Clmtuv});
          }
        }
      }
      /* I have no idea if this will work for the general case, but it works up to d */
      double scale = 0.0;
      for (auto& ituv : tuv)
        scale += ituv.second*ituv.second;

      for (auto& ituv : tuv)
        ituv.second /= scale;

      mtuv.push_back(tuv);
    }
    lmtuv_.push_back(mtuv);
  }
}


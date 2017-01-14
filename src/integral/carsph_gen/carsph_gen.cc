//
// Author : Toru Shiozaki
// Date   : May 2009
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include <boost/lexical_cast.hpp>
#include <cassert>
#include "mpreal.h"

using namespace std;
using namespace mpfr;
using namespace boost;

#define LARGE 32
#define LEND 8

struct Data {
    vector<mpreal> factorial;
    Data() {
      factorial.reserve(LARGE * 2);
      const mpreal one = "1,0";
      factorial.push_back(one);
      for (int i = 1; i != LARGE * 2; ++i) {
        const mpreal imp = i;
        factorial.push_back(factorial.back() * imp);
      }
    };
    ~Data(){};
    const mpreal comb(const int i, const int j) const {
      return factorial[i] / factorial[j] / factorial[i - j];
    };
};


int main() {
  mpfr::mpreal::set_default_prec(10000);
  const mpreal one = "1.0";
  const mpreal two = "2.0";
  const mpreal quarter = "0.25";

  struct Data data;
  vector<mpreal> factorial = data.factorial;

  vector<pair<int, mpreal>> s0(1, make_pair(0, one));
  vector<vector<pair<int, mpreal>>> s1(1, s0);
  vector<vector<vector<pair<int, mpreal>>>> lmtuv(1, s1);

  for (int l = 1; l != LEND; ++l) {

    map<int, int> mapping;
    int cnt = 0;
    for (int z = 0; z <= l; ++z) {
      for (int y = 0; y <= l - z; ++y) {
        const int x = l - y - z;
        if (x < 0) continue;
        const int key = x + y * LARGE + z * LARGE * LARGE;
        mapping.insert(make_pair(key, cnt));
        ++cnt;
      }
    }

    vector<vector<pair<int, mpreal>>> mtuv;
    for (int n = 0; n != 2 * l + 1; ++n) {
      int m = l - (n / 2);
      if (n % 2 == 1) m *= -1;

      const int vm2 = m < 0 ? 1 : 0;
      const int absm = m > 0 ? m : -m;
      vector<pair<int, mpreal>> tuv;

      const mpreal Nlms = one / pow(two, absm) / factorial[l] * sqrt((m == 0 ? one : two) * factorial[l + absm] * factorial[l - absm]);
      const int tmax = floor((l - absm) / 2.0);
      for (int t = 0; t <= tmax; ++t) {
        for (int u = 0; u <= t; ++u) {
          const int vmax = 2 * floor((absm - vm2) / 2.0) + vm2;
          for (int v2 = vm2; v2 <= vmax; v2 += 2) {
            assert((v2 - vm2) % 2 == 0);
            const mpreal Clmtuv = pow(-one, t + ((v2 - vm2) / 2)) * pow(quarter, t)
                                * data.comb(l, t) * data.comb(l - t, absm + t)
                                * data.comb(t, u) * data.comb(absm, v2) * Nlms;
            const int xexp = 2 * t + absm - 2 * u - v2;
            const int yexp = 2 * u + v2;
            const int zexp = l - 2 * t - absm;
            mpreal denom = one;
/*          for (int x = 2; x <= xexp; ++x) denom *= two * x - one;
            for (int y = 2; y <= yexp; ++y) denom *= two * y - one;
            for (int z = 2; z <= zexp; ++z) denom *= two * z - one;
*/
            denom = sqrt(denom);
            map<int, int>::const_iterator current = mapping.find(xexp + yexp * LARGE + zexp * LARGE * LARGE);
            assert(current != mapping.end());
            const mpreal coeff = Clmtuv / denom;
            tuv.push_back(make_pair(current->second, coeff));
          }
        }
      }
      mtuv.push_back(tuv);
    }
    lmtuv.push_back(mtuv);
  }

  ///////////////////////
  // code generator part
  ///////////////////////

// I need to implement S series separately...
  const int cartesian_xyz[20] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66,
                                 78, 91, 105, 120, 136, 153, 171, 190, 210};
  const int spherical_xyz[20] = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21,
                                 23, 25, 27, 29, 31, 33, 35, 37, 39};

  for (int l0 = 0; l0 != LEND; ++l0) {
    const vector<vector<pair<int, mpreal>>> l0mtuv = lmtuv.at(l0);
    assert(l0mtuv.size() == spherical_xyz[l0]);
    const string l0str = lexical_cast<string>(l0);
    for (int l1 = l0; l1 != LEND; ++l1) {

      const string l1str = lexical_cast<string>(l1);
      string filename = "_carsph_" + l1str + l0str + ".cc";
  string code = "\
//\n\
// BAGEL - Brilliantly Advanced General Electronic Structure Library\n\
// Filename: " + filename + "\n\
// Copyright (C) 2009 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// This program is free software: you can redistribute it and/or modify\n\
// it under the terms of the GNU General Public License as published by\n\
// the Free Software Foundation, either version 3 of the License, or\n\
// (at your option) any later version.\n\
//\n\
// This program is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU General Public License\n\
// along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\
//\n\
\n";
    if (l0 == 7 || l1 == 7) {
      code += "\
#ifdef COMPILE_J_ORB\n";
    }
    code += "\
#include <src/integral/carsphlist.h>\n\
#include <algorithm>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n";
  code += "\n";
      const vector<vector<pair<int, mpreal>>> l1mtuv = lmtuv.at(l1);
      assert(l1mtuv.size() == spherical_xyz[l1]);

auto function = [&](const string T, const string classname) {
string code = "\
void " + classname + "::carsph_" + l1str + l0str + "(const int nloop, const " + T + "* source, " + T + "* target) {\n";
      if (l1 > 1) {
        vector<vector<pair<int, mpreal>>>::const_iterator m0;
        vector<vector<pair<int, mpreal>>>::const_iterator m1;
        map<double, string> coeff_map;
        int index = 0;
        for (m1 = l1mtuv.begin(); m1 != l1mtuv.end(); ++m1) {
          for (m0 = l0mtuv.begin(); m0 != l0mtuv.end(); ++m0) {
            for (vector<pair<int, mpreal>>::const_iterator p1 = m1->begin(); p1 != m1->end(); ++p1) {
              for (vector<pair<int, mpreal>>::const_iterator p0 = m0->begin(); p0 != m0->end(); ++p0) {
                const double current = abs(p0->second * p1->second).toDouble();
                if (coeff_map.end() == coeff_map.find(current) && current != 1.0) {
                  const string name = "c" + lexical_cast<string>(index);
                  coeff_map.insert(make_pair(current, name));
                  ++index;
                }
              }
            }
          }
        }
        for (map<double, string>::reverse_iterator iter = coeff_map.rbegin(); iter != coeff_map.rend(); ++iter) {
          code += "  const double " + iter->second + " = " + lexical_cast<string>(iter->first) + ";\n";
        }

        code += "\
  for (int iloop = 0; iloop != nloop; ++iloop, target += " + lexical_cast<string>(l0mtuv.size() * l1mtuv.size())+ ", source += " + lexical_cast<string>(cartesian_xyz[l0] * cartesian_xyz[l1]) + ") {\n\
";
        int cnt = 0;
        for (m1 = l1mtuv.begin(); m1 != l1mtuv.end(); ++m1) {
          for (m0 = l0mtuv.begin(); m0 != l0mtuv.end(); ++m0, ++cnt) { // loops for basis
            code += "\
    target[" + lexical_cast<string>(cnt) + "] = ";
            bool first = true;
            int iformat = 0;
            for (vector<pair<int, mpreal>>::const_iterator p1 = m1->begin(); p1 != m1->end(); ++p1) {
              for (vector<pair<int, mpreal>>::const_iterator p0 = m0->begin(); p0 != m0->end(); ++p0, ++iformat) {
                const int position = cartesian_xyz[l0] * p1->first + p0->first;
                const string posstr = lexical_cast<string>(position);
                const mpreal coeffmp = p0->second * p1->second;
                const double coeff = coeffmp.toDouble();
                const double abscoeff = abs(coeffmp).toDouble();
                const bool abscoeffisone = abscoeff == 1.0;
                map<double, string>::const_iterator found;
                if (!abscoeffisone) {
                  found = coeff_map.find(abscoeff);
                  assert(found != coeff_map.end());
                }
                const string coeffstr = (coeff < 0.0 ? "- " : "") + (!abscoeffisone ? found->second + " * " : "");

                code += (((!first) && coeff >= 0.0) ?  " + " : " ") + coeffstr + "source[" + posstr + "]";
                if (iformat % 3 == 2 && iformat != m0->size() * m1->size() - 1) code += "\n                 ";
                first = false;
              }
            }
            code += ";\n";
          }
        }

      code += "\
  }\n";
      } else {
        code += "\
  copy_n(source, nloop" + (l1 != 0 ? ("*" + lexical_cast<string>(cartesian_xyz[l0] * cartesian_xyz[l1])) : "") + ", target);\n";
      }
    code += "\
}\n\n";
  return code;
};
    code += function("double", "CarSphList");
    code += function("complex<double>", "CCarSphList");

    if (l0 == 7 || l1 == 7) {
      code += "\
#endif\n";
    }

    ofstream ofs;
    ofs.open(filename.c_str());
    ofs << code;
    ofs.close();
    }
  }

  return 0;
}

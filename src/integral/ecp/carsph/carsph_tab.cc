//
// Author : Hai-Anh Le
// Email  : anh@u.northwestern.edu
// Date   : Feb 2014
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
#define LEND 10

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
  const mpreal zero = "0.0";
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

  /////////////////////////////
  //   code generator part   //
  /////////////////////////////

  for (int l1 = 0; l1 != LEND; ++l1) {
    const string l1str = lexical_cast<string>(l1);
    string filename = "_lmtuv_" + l1str + ".cc";
    string code = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: " + filename + "\n\
// Copyright (C) 2009 Toru Shiozaki\n\
//\n\
// Author: Hai-Anh Le <anh@u.northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 3, or (at your option)\n\
// any later version.\n\
//\n\
// The BAGEL package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the BAGEL package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\
\n\
#include <algorithm>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n";
  const int nsph = 2*l1 + 1;
  const int ncar = static_cast<int>(data.comb(l1 + 2, l1).toDouble());
  const string nsphs = lexical_cast<string>(nsph);
  const string ncars = lexical_cast<string>(ncar);
  code += "\n";
      const vector<vector<pair<int, mpreal>>> l1mtuv = lmtuv.at(l1);
      assert(l1mtuv.size() - nsph == 0);

auto function = [&](const string T, const string classname) {

string code = "\
void " + classname + "::lmtuv_" + l1str + "(" + nsphs + ", " + ncars + ") {\n";
code += "constexpr double clmtuv[" + lexical_cast<string>(ncar * nsph) + "] = {";
string indent("   ");
const double tiny = 1.0e-100;

      if (l1 > 1) {
        vector<vector<pair<int, mpreal>>>::const_iterator m1;

        for (m1 = l1mtuv.begin(); m1 != l1mtuv.end(); ++m1) {
          for (vector<pair<int, mpreal>>::const_iterator p1 = m1->begin(); p1 != m1->end(); ++p1) {
            const int position = p1->first;
            const double coeff = (p1->second).toDouble();

            stringstream coeffss;
            coeffss << indent << scientific << setprecision(5) << ((coeff > 0.0 || fabs(coeff) < tiny) ? " " : "")  << setw(10) <<
             (fabs(coeff) < tiny ? 0.0 : coeff);
            code += coeffss.str() + ", ";
          }
          code += ";\n";
        }

      code += "\
  }\n";
      } else {
        code += "\
  copy_n(source, nloop" + (l1 != 0 ? ("*" + lexical_cast<string>(ncar)) : "") + ", target);\n";
      }
    code += "\
}\n\n";
  return code;
};
    code += function("double", "CarSph");

    ofstream ofs;
    ofs.open(filename.c_str());
    ofs << code;
    ofs.close();
    }

  return 0;
}

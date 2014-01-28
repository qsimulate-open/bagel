//
// Author : Hai-Anh Le
// Date   : January 25, 2014
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

constexpr int JMAX = 6;
constexpr int LARGE = 32;

struct Data {
    vector<mpreal> factorial;
    Data() {
      factorial.reserve(LARGE * 2);
      const mpreal one = "1.0";
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

  const int n3j = static_cast<int>((JMAX+1)*(JMAX+1)*(JMAX+1)*(2*JMAX+1)*(2*JMAX+1));
  cout << "n3j = " << n3j << endl;
  vector<double> w3j(n3j);
  struct Data data;
  vector<mpreal> factorial = data.factorial;

  int index = 0;
  int nnzeros = 0;
  const mpreal zero = "0.0";
  for (int j1 = 0; j1 != JMAX + 1; ++j1) {
    for (int j2 = 0; j2 != JMAX + 1; ++j2) {
      for (int j3 = 0; j3 != JMAX + 1; ++j3) {
        for (int ma = 0; ma != 2*JMAX+1; ++ma) {
          for (int mb = 0; mb != 2*JMAX+1; ++mb) {
            const int m1 = ma - JMAX;
            const int m2 = mb - JMAX;
            const int m3 = -(m1 + m2);
            if (j1 + j2 < j3 || j2 + j3 < j1 || j3 + j1 < j2) {
              w3j[index] = 0.0;
            } else {
              if (fabs(m1) > j1 || fabs(m2) > j2 || fabs(m3) > j3) {
              } else {
                if (m3 > j3) {
                  w3j[index] = 0.0;
                } else {
                  w3j[index] = 1.0;
                  mpreal bdel = sqrt(factorial[j1+j2-j3]*factorial[j1-j2+j3]*factorial[-j1+j2+j3]/factorial[j1+j2+j3+1]);
                  mpreal coef = pow(-1, j1-j2-m3);
                  mpreal sqr = sqrt(factorial[j1+m1]*factorial[j1-m1]*factorial[j2+m2]*factorial[j2-m2]*factorial[j3+m3]*factorial[j3-m3]);
                  int kmax = 0;
                  if (kmax < -j3 + j2 - m1) kmax = -j3 + j2 - m1;
                  if (kmax < -j3 + j1 + m2) kmax = -j3 + j1 + m2;
                  int kmin = j1 + j2 - j3;
                  if (kmin > j1 - m1) kmin = j1 - m1;
                  if (kmin > j2 + m2) kmin = j2 + m2;
                  mpreal sumk = zero;
                  if (kmin <= kmax) {
                     for (int k = kmin; k != kmax + 1; ++k) {
                       mpreal a = pow(-1.0, k)/factorial[k]/factorial[j1+j2-j3-k]/factorial[j1-m1-k]/factorial[j2+m2-k];
                       mpreal b = 1.0/factorial[j3-j2+m1+k]/factorial[j3-j1-m2+k];
                       sumk += a*b;
                     }
                  }
                  const mpreal jj = bdel * coef * sqr * sumk;
                  w3j[index] = jj.toDouble();
                  if (w3j[index] != 0.0) ++nnzeros;
                }
              }
            }
//          cout << "index = " << index << "  (" << j1 << ", " << m1 << ") (" <<  j2 << ", " << m2 << ") (" <<  j3 << ", " << m3 << ")" << "  w3j = " << w3j[index] << endl;
            ++index;
          }
        }
      }
    }
  }
  cout << "No. of nonzero terms = " << nnzeros << endl;



  /////////////////////////////////
  // *** code generator part *** //
  /////////////////////////////////

  ofstream ofs;
  const string filename = "wigner3j.h";

  ofs.open(filename.c_str());
  ofs << "\
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
#ifndef __SRC_INTEGRAL_ECP_WIGNER3J_GEN_WIGNER3J_H\n\
#define __SRC_INTEGRAL_ECP_WIGNER3J_GEN_WIGNER3J_H\n\
\n\
#include <algorithm>\n\
#include <cassert>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n\
namespace bagel {\n\
\n\
struct Wigner3j {\n\
  \n\
  public:\n\
  Wigner3j() {}\n\
\n\
  const double lookup_wigner3j(const int j1, const int m1, const int j2, const int m2, const int j3, const int m3) {\n\
  \n";

    ofs << "\
    constexpr double w3j[" << n3j <<"] = {";

    const double tiny = 1.0e-100;
    stringstream listj;
    string indent("   ");
    int cnt = 0;
    for (auto iter = w3j.begin(); iter != w3j.end(); ++iter) {
      listj << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
             (fabs(*iter) < tiny ? 0.0 : *iter);
      if (iter + 1 != w3j.end()) listj << ",";
      if (cnt++ % 7 == 4) listj << "\n";
    }

    ofs << listj.str() << "\
    };" << endl;

    ofs << "\
    \n\
    assert (j1 <= " << JMAX << " && j2 <= " << JMAX << " && j3 <= " << JMAX << "); \n\
    const int m = m1 + m2 + m3; \n\
    if (m != 0) {\n\
      return 0.0;\n\
    } else {\n\
      const int j = j1*" << (JMAX+1)*(JMAX+1)*(2*JMAX+1)*(2*JMAX+1) << " + j2*" << (JMAX +1)*(2*JMAX+1)*(2*JMAX+1) << " + j3*" << (2*JMAX+1)*(2*JMAX+1) << " + (m1+" << JMAX << ")*" << 2*JMAX + 1 << " + m2+" << JMAX << ";\n\
      return w3j[j];\n\
    }\n";

    ofs << "\
    \n\
  }\n\
  \n\
};\n\
}\n\
\n\
#endif\n";

  ofs.close();

  // *** end of code generator part *** //

  return 0;
}

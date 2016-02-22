//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: May 2014
//

#include <src/integral/ecp/sphharmonics.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <memory>
#include <array>
#include <boost/lexical_cast.hpp>

using namespace bagel;
using namespace std;
using namespace boost;

constexpr int LMAX = 10;

int main() {

  vector<vector<vector<double>>> usp;

  for (int l = 0; l <= LMAX; ++l) {
    vector<vector<double>> lusp;
    const size_t size = (l+1)*(l+2)/2;

    for (int mu = 0; mu <= 2*l; ++mu) {
      const int m = mu - l;
      vector<double> lmusp;

      array<int, 2> lm = {l, m};
      shared_ptr<SphHarmonics> sphusp = make_shared<SphHarmonics>(lm);
      for (int lz = 0; lz <= l; ++lz) {
        for (int ly = 0; ly <= l - lz; ++ly) {
          const int lx = l - lz - ly;
          const double coeff = sphusp->sph_to_USP(lx, ly);
          lmusp.push_back(coeff);
        }
      }

      assert(lmusp.size() == size);

      lusp.push_back(lmusp);
    }

    usp.push_back(lusp);
  }

  ////////////////////////
  // code generator part /
  ////////////////////////

  for (int l = 0; l <= LMAX; ++l) {
    ofstream ofs;
    const string lstr = lexical_cast<string>(l);
    string filename = "_sphusp_" + lstr + ".cc";
    const int size = (l + 1) * (l + 2) * (2*l + 1) / 2;

  ofs.open(filename.c_str());
  ofs << "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: " + filename + "\n\
// Copyright (C) 2014 Toru Shiozaki\n\
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
#include <cassert>\n\
#include <src/integral/ecp/sphusplist.h>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n";

ofs << "\
vector<double> SphUSPList::sphusp_" << lstr << "(const int m) {\n";
  ofs << "\n\
  vector<double> c;\n\
  constexpr double coeff[" << size << "] = {";

    const vector<vector<double>> lusp = usp.at(l);
    stringstream listc;
#if 1
    for (int mu = 0; mu <= 2*l; ++mu) {
      const vector<double> lmusp = lusp.at(mu);

      const double tiny = 1.0e-100;
      int cnt = 0;

      string indent("  ");
      for (auto iter = lmusp.begin(); iter != lmusp.end(); ++iter) {
        listc << indent << scientific << setprecision(15) << ((*iter > 0.0 || fabs(*iter) < tiny) ? " " : "")  << setw(20) <<
                 (fabs(*iter) < tiny ? 0.0 : *iter);
        if (iter + 1 != lmusp.end() || mu != 2*l) listc << ",";
        if (cnt++ % 5 == 4) listc << "\n";
      }

    }
#endif

  ofs << listc.str() << "};\n\
  \n";

  ofs << "\
  assert(abs(m) <= " << l << ");\n\
  const int size_c = (" << l << " + 1) * (" << l << " + 2) / 2;\n\
  const int mu = m + " << l << ";\n\
  const int i0 = mu * size_c;\n\
  for (int i = i0; i != i0 + size_c; ++i) c.push_back(coeff[i]);\n\
  return c;\n\
  \n";

ofs << "}\n";

ofs.close();

  }

  return 0;

}

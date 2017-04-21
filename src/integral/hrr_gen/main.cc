//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: main.cc
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


#include <string>
#include <fstream>
#include "../../util/constants.h"
#include "angular_index.h"
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <list>
#include <vector>
#include <algorithm>

using namespace std;
using namespace boost;
using namespace bagel;

int main() {

  for (int n = 0; n != ANG_VRR_END; ++n) {
    for (int a = 0; a <= n; ++a) {
      const int b = n - a;
      if (a < b || a >= ANG_HRR_END) continue;
      if (b == 0) continue;
      string astr;
      if (a < 10) astr = lexical_cast<string>(a);
      if (a == 10) astr = lexical_cast<string>("a");
      if (a == 11) astr = lexical_cast<string>("b");
      if (a == 12) astr = lexical_cast<string>("c");
      if (a == 13) astr = lexical_cast<string>("d");
      if (a == 14) astr = lexical_cast<string>("e");
      string nstr;
      if (n < 10) nstr = lexical_cast<string>(n);
      if (n == 10) nstr = lexical_cast<string>("a");
      if (n == 11) nstr = lexical_cast<string>("b");
      if (n == 12) nstr = lexical_cast<string>("c");
      if (n == 13) nstr = lexical_cast<string>("d");
      if (n == 14) nstr = lexical_cast<string>("e");

      string n0 = nstr + lexical_cast<string>(0);
      string ab = astr + lexical_cast<string>(b);
      string sa = astr;
      string out;
      const double cartesian_xyz[20] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66,
                                       78, 91, 105, 120, 136, 153, 171, 190, 210};

      int size = 0;
      for (int i = max(a, b); i <= a + b; ++i) size += cartesian_xyz[i];

      int mapping[(n + 1) * (n + 1) * (n + 1)];
      int cnt_local = 0;
      for (int i = max(a,b); i <= a + b; ++i) {
        for (int iz = 0; iz <= i; ++iz) {
          for (int iy = 0; iy <= i - iz; ++iy) {
            const int ix = i - iy - iz;
            if (ix >= 0) {
              mapping[ix + (n + 1) * (iy + (n + 1) * iz)] = cnt_local;
             ++cnt_local;
            }
          }
        }
      }

      string ssize = lexical_cast<string>(size);

      size = cartesian_xyz[a] * cartesian_xyz[b];
      string osize = lexical_cast<string>(size);
      size = cartesian_xyz[a];
      string asize = lexical_cast<string>(size);

      string filename = "_hrr_" + n0 + "_" + ab + ".cc";
      out +="\
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
      if (a == 7 || b == 7) {
        out += "\
#ifdef COMPILE_J_ORB\n";
      }
      out +="\
#include <src/integral/hrrlist.h>\n\
#include <array>\n\
#include <algorithm>\n\
\n\
using namespace std;\n\
using namespace bagel;\n\
\n";
      auto function = [&](const string T, const string classname) {
        string out = "\
void " + classname + "::perform_HRR_" + n0 + "_" + ab + "(const int nloop, const " + T + "* data_start, const array<double,3>& AB, " + T + "* data_out) {\n\
  for (int c = 0; c != nloop; ++c) {\n\
    auto current_data = &data_start[c * " + ssize + "];\n\
    auto current_out = &data_out[c * " + osize + "];\n";

    int cnt = 0;
    string text0, text1, text2, text3;


    int number = 0;
    for (int jz = 0; jz <= a; ++jz) {
      const int rjz = (n + 1) * jz;
      for (int jy = 0; jy <= a - jz; ++jy) {
        const int jx = a - jy - jz;
        if (jx >= 0) {
          double pl = jx + (n + 1) * (rjz + jy);
          Angular_Index base(jx, jy, jz);
      out += "\
   {\n\
     //current index a: " + base.show() + "\n";
      string out1, out2;
      for (int m = 0; m <= (n - a); ++m) {
        for (int iz = 0; iz <= m; ++iz) {
          for (int iy = 0; iy <= m - iz; ++iy) {

            const int ix = m - iz - iy;
            if (ix >= 0) {
              Angular_Index ai(ix, iy, iz);
              const string offset = lexical_cast<string>(pl + ix + (n + 1) * (iy + (n + 1) * iz));
              const int intoffset =  pl + ix + (n + 1) * (iy + (n + 1) * iz);
              const int mapinto = mapping[intoffset];
//            out1 += "\
//    const int ja" + ai.show() + "_0 = mapping[" + offset + "];\n";
              out2 += "\
      const auto a" + ai.show() + "_0 = current_data[" + lexical_cast<string>(mapinto) + "];\n";
//    const double a" + ai.show() + "_0 = current_data[ja" + ai.show() + "_0];\n";
            }
          }
        }
      }
      out += out2 + "\n";
//    out += out1 + "\n" + out2;

      list<Angular_Pair> needed;
      vector<string> blocks;

      for (int m = b; m >= 1; --m) {
        typedef list<string> AP;
        vector<AP> mylist;
        mylist.resize(ANG_HRR_END);
        for (int iz = 0; iz <= m; ++iz) {
          for (int iy = 0; iy <= m - iz; ++iy) {
            const int ix = m - iz - iy;
            if (ix >= 0) {
              Angular_Index bi(ix, iy, iz);
              for (int l = n - m - a; l >= 0; --l) {
                for (int jz = 0; jz <= l; ++jz) {
                  for (int jy = 0; jy <= l - jz; ++jy) {
                    const int jx = l - jz - jy;
                    if (jx >= 0) {
                      Angular_Index ai(jx, jy, jz);
                      Angular_Pair indexpair(make_pair(ai, bi));
                      if (m == b || find(needed.begin(), needed.end(), indexpair) != needed.end()) {
                        tuple<Angular_Pair, Angular_Pair, int> next = indexpair.hrr_formula();
                        needed.push_back(get<0>(next));
                        needed.push_back(get<1>(next));
                        const int ixyz = jx + jy + jz;
                        string number_str = lexical_cast<string>(number);
                        if (m != b) {
                          const string sprim = "\
      const auto " + indexpair.show() + " = " + get<0>(next).show() +
                        " + AB[" + lexical_cast<string>(get<2>(next)) + "] * " + get<1>(next).show() + ";\n";
                          mylist[ixyz].push_back(sprim);
                        } else {
                          const string sprim = "\
      current_out[" + lexical_cast<string>(cnt) + "] = " + get<0>(next).show() +
                        " + AB[" + lexical_cast<string>(get<2>(next)) + "] * " + get<1>(next).show() + "; // " + indexpair.show() + "\n";
                          mylist[ixyz].push_back(sprim);
                          ++cnt;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        string outblock;
        for (vector<AP>::iterator iter = mylist.begin(); iter != mylist.end(); ++iter) {
          bool goinside = false;
          for (AP::iterator aiter = iter->begin(); aiter != iter->end(); ++aiter) {
            goinside = true;
            outblock += *aiter;
          }
          if (goinside) outblock += "\n";
        }
        blocks.push_back(outblock);
      }
      for (vector<string>::reverse_iterator riter = blocks.rbegin(); riter != blocks.rend(); ++riter)
        out += *riter;
      out += "\
    }\n";

    }}}
      out +="\
  }\n\
}\n";
        return out;
      }; // end of function
      out += function("double", "HRRList");
      out += "\n\n";
      out += function("complex<double>", "CHRRList");
      if (a == 7 || b == 7) {
        out += "\
#endif\n";
      }
      ofstream ofs;
      ofs.open(filename.c_str());
      ofs << out << endl;
      ofs.close();
    }
  }
  return 0;
}

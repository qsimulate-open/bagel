//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <string>
#include <iomanip>
#include <sstream>

#include <fstream>
#define SPHERICAL 1

constexpr int ANG_HRR_END = 8;

using namespace std;

template<typename T>
string String(T a) {
  stringstream ss;
  ss << a;
  return ss.str();
};

int main() {

  ofstream ofs;
#if SPHERICAL
  string filename = "csort_sph.cc";
#else
  string filename = "csort.cc";
#endif
  ofs.open(filename.c_str());

  stringstream out;
  out << "//" << endl;
  out << "// BAGEL - Brilliantly Advanced General Electronic Structure Library" << endl;
  out << "// Filename: " << filename << endl;
  out << "// Copyright (C) 2009 Toru Shiozaki" << endl;
  out << "//" << endl;
  out << "// Author: Toru Shiozaki <shiozaki@northwestern.edu>" << endl;
  out << "// Maintainer: Shiozaki group" << endl;
  out << "//" << endl;
  out << "// This file is part of the BAGEL package." << endl;
  out << "//" << endl;
  out << "// This program is free software: you can redistribute it and/or modify" << endl;
  out << "// it under the terms of the GNU General Public License as published by" << endl,
  out << "// the Free Software Foundation, either version 3 of the License, or" << endl;
  out << "// (at your option) any later version." << endl;
  out << "//" << endl;
  out << "// This program is distributed in the hope that it will be useful," << endl;
  out << "// but WITHOUT ANY WARRANTY; without even the implied warranty of" << endl;
  out << "// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the" << endl;
  out << "// GNU General Public License for more details." << endl;
  out << "//" << endl;
  out << "// You should have received a copy of the GNU General Public License" << endl;
  out << "// along with this program.  If not, see <http://www.gnu.org/licenses/>." << endl;
  out << "//" << endl;
  out << endl;
  out << endl;
  out << "#include <algorithm>" << endl;
  out << "#include <src/integral/sortlist.h>" << endl;
  out << endl;
  out << "using namespace std;" << endl;
  out << "using namespace bagel;" << endl;
  out << endl;

  for (int x2 = 0; x2 != ANG_HRR_END; ++x2) {
  for (int x3 = 0; x3 != x2 + 1; ++x3) {
#if SPHERICAL
    const double xyz[10] = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19};
#else
    const double xyz[10] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55};
#endif
    const int x2end = xyz[x2];
    const int x3end = xyz[x3];
    const string label = String<int>(x3) + String<int>(x2);

    if (x2 == 7 && x3 == 0) {
      out << "\
#ifdef COMPILE_J_ORB\n";
    }

#if SPHERICAL
    out << "\
void CSortList::sort_indices_" << label << "_sph(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {\n";
#else
    out << "\
void CSortList::sort_indices_" << label << "(complex<double>* target, const complex<double>* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {\n";
#endif
    out << "\
  const int innerloopsize = c2end * c3end * " << x2end*x3end << ";\n\
  if (!swap23) {\n\
    int offset = 0;\n\
    const int cont2csize = " << x2end << " * c2end;\n\
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {\n\
      complex<double>* current_target = &target[offset];\n\
      const complex<double>* current_source = &source[offset];\n\
\n\
      for (int c2 = 0; c2 != c2end; ++c2) {\n\
        const int c2x2end = c2 * " << x2end << ";\n\
        for (int c3 = 0; c3 != c3end; ++c3) {\n\
          const int soffset = " << x2end*x3end << " * (c3 + c3end * c2);\n\
          const int toffset = " << x3end << " * c3 * cont2csize + c2x2end;" << endl;
    for (int x2 = 0; x2 != x2end; ++x2) {
      for (int x3 = 0; x3 != x3end; ++x3) {
        out << "\
          current_target[toffset + " << x3 << " * cont2csize + " << x2 << "] = current_source[soffset + " << (x3 + x3end * (x2)) << "];" << endl;
      }
    }
    out << "\
        }\n\
      }\n\
\n\
    }\n\
  } else {\n\
    int offset = 0;\n\
    const int cont3csize = " << x3end << " * c3end;\n\
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {\n\
      complex<double>* current_target = &target[offset];\n\
      const complex<double>* current_source = &source[offset];\n\
\n\
      for (int c2 = 0; c2 != c2end; ++c2) {\n\
        for (int c3 = 0; c3 != c3end; ++c3) {\n\
          const int c3x3end = c3 * " << x3end << ";\n";
          out << "\
          const int soffset = " << x2end*x3end << " * (c3 + c3end * c2);\n\
          const int toffset = " << x2end << " * c2 * cont3csize + c3x3end;" << endl;
          for (int x2 = 0; x2 != x2end; ++x2) {
            out << "\
          copy_n(current_source+soffset+" << setw(3) << (x3end * (x2)) << ", " << setw(3) << x3end << ", current_target+toffset+" << setw(2) << x2 << "*cont3csize);" << endl;
          }
    out << "\
        }\n\
      }\n\
\n\
    }\n\
  }\n\
\n\
}" << endl << endl << endl;

}}

  if (ANG_HRR_END > 7) {
    out << "\
#endif\n";
  }

  ofs << out.str() << endl;
  ofs.close();
  return 0;
}

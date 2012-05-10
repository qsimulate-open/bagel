//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <string>
#include <sstream>
#include "../macros.h"
  
#include <fstream>
//#define SPHERICAL 1

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
  string filename = "sort_sph.cc";
#else
  string filename = "sort.cc";
#endif
  ofs.open(filename.c_str());

  stringstream out;
    out << "\
//\n\
// Author: Toru Shiozaki\n\
// Machine Generated Code\n\
//\n\
\n\
#include \"sortlist.h\"\n\
#include <algorithm>\n\
\n\
using namespace std;\n\
\n";

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

#if SPHERICAL
    out << "\
void SortList::sort_indices_" << label << "_sph(double* target, const double* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {\n";
#else
    out << "\
void SortList::sort_indices_" << label << "(double* target, const double* source, const int c3end, const int c2end, const int loopsize, const bool swap23) {\n";
#endif
    out << "\
  const int innerloopsize = c2end * c3end * " << x2end*x3end << ";\n\
  if (!swap23) {\n\
    int offset = 0;\n\
    const int cont2csize = " << x2end << " * c2end; \n\
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {\n\
      double* current_target = &target[offset];\n\
      const double* current_source = &source[offset];\n\
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
      double* current_target = &target[offset];\n\
      const double* current_source = &source[offset];\n\
\n\
      for (int c2 = 0; c2 != c2end; ++c2) {\n\
        for (int c3 = 0; c3 != c3end; ++c3) {\n\
          const int c3x3end = c3 * " << x3end << ";\n";
          out << "\
          const int soffset = " << x2end*x3end << " * (c3 + c3end * c2);\n\
          const int toffset = " << x2end << " * c2 * cont3csize + c3x3end;" << endl; 
          for (int x2 = 0; x2 != x2end; ++x2) {
            out << "\
          copy(current_source+soffset+" << (x3end * (x2)) << ", current_source+soffset+" << (x3end * (x2+1)) << ", current_target+toffset+" << x2 << "*cont3csize);" << endl; 
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
  ofs << out.str() << endl;
  ofs.close();
  return 0;
}

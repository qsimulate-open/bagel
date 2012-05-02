//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <iostream>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include "../../rysint/macros.h"
#include "vrr.h"

using namespace std;
using namespace boost;

int main() {

  vector<string> declist;

  for (int i = 0; i != ANG_VRR_END; ++i) {
    for (int j = 0; j != ANG_VRR_END; ++j) {
#if 0
      if (i == 0 && j == 0) continue;
      if (i == 1 && j == 0) continue;
      if (i == 0 && j == 1) continue;
//sum 2
      if (i == 2 && j == 0) continue;
      if (i == 0 && j == 2) continue;
      if (i == 1 && j == 1) continue;
//sum 3
      if (i == 3 && j == 0) continue;
      if (i == 0 && j == 3) continue;
      if (i == 2 && j == 1) continue;
      if (i == 1 && j == 2) continue;
//sum 4
//    if (i == 4 && j == 0) continue;
//    if (i == 0 && j == 4) continue;
      if (i == 2 && j == 2) continue;
#endif
      VRR vrrij(i, j);
      string istr, jstr;
      if (i < 10) istr = lexical_cast<string>(i);
      else if (i == 10) istr = lexical_cast<string>("a");
      else if (i == 11) istr = lexical_cast<string>("b");
      else if (i == 12) istr = lexical_cast<string>("c");
      else if (i == 13) istr = lexical_cast<string>("d");
      else if (i == 14) istr = lexical_cast<string>("e");
      else if (i == 15) istr = lexical_cast<string>("f");
      if (j < 10) jstr = lexical_cast<string>(j);
      else if (j == 10) jstr = lexical_cast<string>("a");
      else if (j == 11) jstr = lexical_cast<string>("b");
      else if (j == 12) jstr = lexical_cast<string>("c");
      else if (j == 13) jstr = lexical_cast<string>("d");
      else if (j == 14) jstr = lexical_cast<string>("e");
      else if (j == 15) jstr = lexical_cast<string>("c");
      const string filename = "_gvrr_" + istr + "0" + jstr + "0"; 
      declist.push_back(vrrij.dump(filename));
//    cout << "    static void " << filename << "(double*, const double*, const double*, const double*, const double*, const double*);" << endl;
      cout << filename << ".cc "; 
    }
  }
  cout << endl;

  return 0;
}


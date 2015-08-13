//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: August 2015
//

#include <iomanip>
#include "localexps.h"

using namespace std;

int main() {

  vector<complex<double>> mlm;
  const int ws = 1;
  const int lmax = 5;
  const int limit = 10;
  const double thresh = 1e-14;
  LocalExps test(ws, lmax, limit, thresh);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l) {
    for (int mm = 0; mm <= 2*l+1; ++mm, ++cnt) {
      const int m = mm - l;
      const double mlm = test.mlm_real(cnt);
      if (abs(mlm) > 1e-10)
        cout << "l = " << l << "  m = " << m << "  mlm = " << setw(20) << setprecision(14) << mlm << endl;
    }
  }

  return 0;
}

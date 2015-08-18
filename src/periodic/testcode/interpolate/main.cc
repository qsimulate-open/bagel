//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: August 2015
//

#include "localexps.h"
#include <src/util/parallel/resources.h>
#include <src/util/parallel/mpi_interface.h>

using namespace std;

namespace bagel{
  Resources* resources__;
  MPI_Interface* mpi__;
}

int main() {
  mpi__ = new MPI_Interface();
  resources__ = new Resources(8);

  vector<complex<double>> mlm;
  const int ws = 1;
  const int lmax = 5;
  const int limit = 2;
  LocalExps test(ws, lmax, limit);
  int cnt = 0;
  for (int l = 0; l <= lmax; ++l) {
    for (int mm = 0; mm <= 2*l+1; ++mm, ++cnt) {
      const int m = mm - l;
      const double mlm = test.mlm_real(cnt);
      if (abs(mlm) > 1e-10)
        cout << "l = " << l << "  m = " << m << "  mlm = " << setw(20) << setprecision(14) << mlm << endl;
    }
  }

  delete resources__;
  delete mpi__;

  return 0;
}

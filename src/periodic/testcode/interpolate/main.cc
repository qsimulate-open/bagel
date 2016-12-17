//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: August 2015
//

#include "localexps.h"
#include <src/util/parallel/resources.h>
#include <src/util/parallel/mpi_interface.h>

using namespace std;
using namespace test;

namespace bagel{
  Resources* resources__;
  MPI_Interface* mpi__;
}

int main() {
  mpi__ = new MPI_Interface();
  resources__ = new Resources(8);

  vector<complex<double>> mlm;
  const int ws = 1;
  const int lmax = 20;
  const int limit = 10;


  LocalExps test(ws, lmax, limit);
  for (int l = 0; l <= lmax; ++l) {
    for (int m = 0; m <= l; ++m) { // Mlm = -Ml-m
      const int imul = l * l + m + l;
      const double mlm = test.mlm_real(imul);
      const double tmp = test.mlm_imag(imul);
      if (abs(mlm) > 1e-8)
        cout << "l = " << l << "  m = " << m << "  mlm = " << setw(20) << setprecision(14) << mlm << "     " << tmp << endl;
    }
  }

  delete resources__;
  delete mpi__;

  return 0;
}

//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <src/osint/f77.h>
#include <src/osint/osint.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <src/stackmem.h>

using namespace std;

extern StackMem* stack;

void OSInt::perform_contraction(const int asize, const double* prim, const int pdim0, const int pdim1, double* cont, 
                                       const vector<vector<double> >& coeff0, const vector<pair<int, int> >& ranges0, const int cdim0, 
                                       const vector<vector<double> >& coeff1, const vector<pair<int, int> >& ranges1, const int cdim1) {
  // transformation of index1
  const int unit = 1;
  const double zero = 0.0;
  const int zeroint = 0;
  const int worksize = pdim1 * asize;
  double* work = stack->get(worksize);

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = ranges0[i].first;
    const int end0   = ranges0[i].second;
    fill(work, work + worksize, zero);
    for (int j = begin0; j != end0; ++j) 
      daxpy_(&worksize, &coeff0[i][j], &prim[j * worksize], &unit, work, &unit); 

    for (int k = 0; k != cdim1; ++k, cont += asize) {
      const int begin1 = ranges1[k].first;
      const int end1   = ranges1[k].second;
      for (int j = begin1; j != end1; ++j) {
        daxpy_(&asize, &coeff1[k][j], &work[j * asize], &unit, cont, &unit); 
      }
    }
  }
  stack->release(worksize);
}



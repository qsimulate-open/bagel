//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/scf/overlap.h>
#include <src/osint/overlapbatch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>

typedef boost::shared_ptr<Geometry> RefGeometry;
typedef boost::shared_ptr<Shell> RefShell;

using namespace std;

Overlap::Overlap(const RefGeometry gm) : Matrix1e(gm) {

  init();
  symmetrize();

}


Overlap::~Overlap() {

}


void Overlap::computebatch(const vector<RefShell>& input, const int offsetb0, const int offsetb1, const int ndim) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis(); 
  const int dimb0 = input[1]->nbasis(); 
  OverlapBatch overlap(input);
  overlap.compute();
  const double* odata = overlap.data();

  int cnt = 0;
  for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
      data_[i * ndim + j] = odata[cnt];
    }
  }
}



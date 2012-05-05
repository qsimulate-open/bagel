//
// Author : Toru Shiozaki
// Date   : May 2012
//

#include <src/grad/gnaibatch.h>

using namespace std;

void GNAIBatch::set_exponents() {
  exponents_ = unique_ptr<double[]>(new double[primsize_*2]);
  double* tmp = exponents_.get();
  for (auto i0 = basisinfo_[0]->exponents().begin(); i0 != basisinfo_[0]->exponents().end(); ++i0) {
    for (auto i1 = basisinfo_[1]->exponents().begin(); i1 != basisinfo_[1]->exponents().end(); ++i1, tmp+=2) {
      tmp[0] = *i0;
      tmp[1] = *i1;
    }
  }
}

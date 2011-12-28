//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include <src/slater/svrrlist.h>

// returns double array of length 6
void SVRRList::_svrr_0020(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;

  data_[2] = D00_[0];
  data_[3] = D00_[1];

  data_[4] = D00_[0] * data_[2] + B01_[0];
  data_[5] = D00_[1] * data_[3] + B01_[1];
}


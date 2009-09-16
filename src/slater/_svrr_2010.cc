//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include "svrrlist.h"

// returns double array of length 18
void SVRRList::_svrr_2010(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;

  data_[3] = C00_[0];
  data_[4] = C00_[1];
  data_[5] = C00_[2];

  data_[6] = C00_[0] * data_[3] + B10_[0];
  data_[7] = C00_[1] * data_[4] + B10_[1];
  data_[8] = C00_[2] * data_[5] + B10_[2];

  data_[9] = D00_[0];
  data_[10] = D00_[1];
  data_[11] = D00_[2];

  data_[12] = C00_[0] * data_[9] + B00_[0];
  data_[13] = C00_[1] * data_[10] + B00_[1];
  data_[14] = C00_[2] * data_[11] + B00_[2];

  data_[15] = C00_[0] * data_[12] + B10_[0] * data_[9] + B00_[0] * data_[3];
  data_[16] = C00_[1] * data_[13] + B10_[1] * data_[10] + B00_[1] * data_[4];
  data_[17] = C00_[2] * data_[14] + B10_[2] * data_[11] + B00_[2] * data_[5];
}


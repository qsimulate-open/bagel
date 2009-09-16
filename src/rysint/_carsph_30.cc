//
// Author : Toru Shiozaki
// Date   : May 2009
//

#include "carsphlist.h"
#include "f77.h"


void CarSphList::carsph_30(const int nloop, const double* source, double* target) {
  const double c3 = 3.872983346207417;
  const double c4 = 2.4494897427831779;
  const double c1 = 2.3717082451262845;
  const double c2 = 1.9364916731037085;
  const double c6 = 1.5;
  const double c0 = 0.79056941504209488;
  const double c5 = 0.61237243569579447;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 7, source += 10) {
    target[0] =  c0 * source[0] - c1 * source[2];
    target[1] =  c1 * source[1] - c0 * source[3];
    target[2] =  c2 * source[4] - c2 * source[6];
    target[3] =  c3 * source[5];
    target[4] =  c4 * source[7] - c5 * source[0] - c5 * source[2];
    target[5] =  c4 * source[8] - c5 * source[1] - c5 * source[3];
    target[6] =  source[9] - c6 * source[4] - c6 * source[6];
  }
}


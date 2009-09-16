//
// Author : Toru Shiozaki
// Date   : May 2009
//

#include "carsphlist.h"
#include "f77.h"


void CarSphList::carsph_22(const int nloop, const double* source, double* target) {
  const double c4 = 3;
  const double c5 = 1.7320508075688772;
  const double c1 = 1.5;
  const double c2 = 0.8660254037844386;
  const double c0 = 0.75;
  const double c6 = 0.5;
  const double c3 = 0.4330127018922193;
  const double c7 = 0.25;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 25, source += 36) {
    target[0] =  c0 * source[0] - c0 * source[2] - c0 * source[12]
                  + c0 * source[14];
    target[1] =  c1 * source[1] - c1 * source[13];
    target[2] =  c1 * source[3] - c1 * source[15];
    target[3] =  c1 * source[4] - c1 * source[16];
    target[4] =  c2 * source[5] - c3 * source[0] - c3 * source[2]
                  - c2 * source[17] + c3 * source[12] + c3 * source[14];
    target[5] =  c1 * source[6] - c1 * source[8];
    target[6] =  c4 * source[7];
    target[7] =  c4 * source[9];
    target[8] =  c4 * source[10];
    target[9] =  c5 * source[11] - c2 * source[6] - c2 * source[8];
    target[10] =  c1 * source[18] - c1 * source[20];
    target[11] =  c4 * source[19];
    target[12] =  c4 * source[21];
    target[13] =  c4 * source[22];
    target[14] =  c5 * source[23] - c2 * source[18] - c2 * source[20];
    target[15] =  c1 * source[24] - c1 * source[26];
    target[16] =  c4 * source[25];
    target[17] =  c4 * source[27];
    target[18] =  c4 * source[28];
    target[19] =  c5 * source[29] - c2 * source[24] - c2 * source[26];
    target[20] =  c2 * source[30] - c2 * source[32] - c3 * source[0]
                  + c3 * source[2] - c3 * source[12] + c3 * source[14];
    target[21] =  c5 * source[31] - c2 * source[1] - c2 * source[13];
    target[22] =  c5 * source[33] - c2 * source[3] - c2 * source[15];
    target[23] =  c5 * source[34] - c2 * source[4] - c2 * source[16];
    target[24] =  source[35] - c6 * source[30] - c6 * source[32]
                  - c6 * source[5] + c7 * source[0] + c7 * source[2]
                  - c6 * source[17] + c7 * source[12] + c7 * source[14];
  }
}


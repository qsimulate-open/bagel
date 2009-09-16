//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <src/rysint/eribatch.h>

using namespace std;

void ERIBatch::sort_indices(double* target, const double* source, const int x3end, const int x2end,
                            const int c3end, const int c2end, const int loopsize, const bool swap23) {
  // from cont2{ cont3{ c{ d{ } } } }
  // to cont3{ d{ cont2{ c{ } } } }     if !swap23_ 
  // to cont2{ c{ cont3{ d{ } } } }     if swap23_ 
  // called by the following format: sort_indices(bkup_, data_, d, c, cont3size_, cont2size_, swap23_);

  const int innerloopsize = c2end * c3end * x2end * x3end;
  if (!swap23) {
    int offset = 0;
    const int cont2csize = x2end * c2end; 
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      double* current_target = &target[offset];

      for (int c2 = 0; c2 != c2end; ++c2) { 
        for (int c3 = 0; c3 != c3end; ++c3) { 
          int i2 = x2end * c2;
          for (int x2 = 0; x2 != x2end; ++x2, ++i2) {
            int ci3 = x3end * c3 * cont2csize;
            for (int x3 = 0; x3 != x3end; ++x3, ci3 += cont2csize, ++source) {
              current_target[ci3 + i2] = *source; 
            }
          }
        }
      }

    }
  } else {
    int offset = 0;
    const int cont3csize = x3end * c3end; 
    for (int i = 0; i != loopsize; ++i, offset += innerloopsize) {
      double* current_target = &target[offset];

      for (int c2 = 0; c2 != c2end; ++c2) { 
        for (int c3 = 0; c3 != c3end; ++c3) { 
          int ci2 = x2end * c2 * cont3csize;
          for (int x2 = 0; x2 != x2end; ++x2, ci2 += cont3csize) {
            int i3 = x3end * c3;
            for (int x3 = 0; x3 != x3end; ++x3, ++i3, ++source) {
              current_target[i3 + ci2] = *source; 
            }
          }
        }
      }

    }
  }

}

//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <src/osint/overlapbatch.h>
#include <iostream>
#include <iomanip>
#include <src/stackmem.h>

using namespace std;

extern StackMem* stack;

// private functions
void OverlapBatch::perform_VRR(double* intermediate) {

  const int worksize = amax1_;
  double* workx = stack->get(worksize);
  double* worky = stack->get(worksize);
  double* workz = stack->get(worksize);

  // Perform VRR
  for (int ii = 0; ii != prim0_ * prim1_; ++ii) {
    const int offset_ii = ii * asize_;
    double* current_data = &intermediate[offset_ii];

    /// Sx(0 : i + j, 0) etc will be made here 
    workx[0] = coeffsx_[ii]; 
    worky[0] = coeffsy_[ii];
    workz[0] = coeffsz_[ii];
    if (ang0_ + ang1_ > 0) {
      workx[1] = (p_[ii * 3    ] - basisinfo_[0]->position(0)) * workx[0];
      worky[1] = (p_[ii * 3 + 1] - basisinfo_[0]->position(1)) * worky[0];
      workz[1] = (p_[ii * 3 + 2] - basisinfo_[0]->position(2)) * workz[0];
      for (int i = 2; i != amax1_; ++i) {
        workx[i] = (p_[ii * 3    ] - basisinfo_[0]->position(0)) * workx[i - 1] + 0.5 * (i - 1) / xp_[ii] * workx[i - 2]; 
        worky[i] = (p_[ii * 3 + 1] - basisinfo_[0]->position(1)) * worky[i - 1] + 0.5 * (i - 1) / xp_[ii] * worky[i - 2]; 
        workz[i] = (p_[ii * 3 + 2] - basisinfo_[0]->position(2)) * workz[i - 1] + 0.5 * (i - 1) / xp_[ii] * workz[i - 2]; 
      } 
    }

    /// assembly process

    for (int iz = 0; iz <= amax_; ++iz) { 
      for (int iy = 0; iy <= amax_ - iz; ++iy) { 
        const double iyiz = workz[iz] * worky[iy]; 
        for (int ix = max(0, amin_ - iy - iz); ix <= amax_ - iy - iz; ++ix) { 
          int pos = amapping_[ix + amax1_ * (iy + amax1_ * iz)];
          current_data[pos] = workx[ix] * iyiz;
        }
      }
    }

  } // end of primsize loop

  stack->release(worksize*3);
}


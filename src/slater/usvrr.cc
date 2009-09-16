//
// Author: Toru Shiozaki
// Machine generated code
//

#include <src/slater/slaterbatch.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cstring>

using namespace std;

/*
void SlaterBatch::perform_USVRR3() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 3 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[3];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 3, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[3];
    double wt[3];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 3, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 3, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}
*/


void SlaterBatch::perform_USVRR4() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 4 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[4];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 4, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[4];
    double wt[4];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 4, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 4, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR5() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 5 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[5];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 5, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[5];
    double wt[5];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 5, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 5, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR6() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 6 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[6];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 6, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[6];
    double wt[6];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    wt[5] = weights_[offset + 5] * roots_[offset + 5];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    womt[5] = weights_[offset + 5]  - wt[5];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 6, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 6, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            iyiz[5] = worky[offsety + 5] * workz[offsetz + 5]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data[ijposition] += iyiz[5] * workx[offsetx + 5]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
                current_data_y[ijposition] += iyiz[5] * workx_y[offsetx + 5]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR7() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 7 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[7];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 7, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[7];
    double wt[7];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    wt[5] = weights_[offset + 5] * roots_[offset + 5];
    wt[6] = weights_[offset + 6] * roots_[offset + 6];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    womt[5] = weights_[offset + 5]  - wt[5];
    womt[6] = weights_[offset + 6]  - wt[6];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 7, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 7, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            iyiz[5] = worky[offsety + 5] * workz[offsetz + 5]; 
            iyiz[6] = worky[offsety + 6] * workz[offsetz + 6]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data[ijposition] += iyiz[5] * workx[offsetx + 5]; 
                current_data[ijposition] += iyiz[6] * workx[offsetx + 6]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
                current_data_y[ijposition] += iyiz[5] * workx_y[offsetx + 5]; 
                current_data_y[ijposition] += iyiz[6] * workx_y[offsetx + 6]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR8() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 8 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[8];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 8, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[8];
    double wt[8];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    wt[5] = weights_[offset + 5] * roots_[offset + 5];
    wt[6] = weights_[offset + 6] * roots_[offset + 6];
    wt[7] = weights_[offset + 7] * roots_[offset + 7];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    womt[5] = weights_[offset + 5]  - wt[5];
    womt[6] = weights_[offset + 6]  - wt[6];
    womt[7] = weights_[offset + 7]  - wt[7];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 8, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 8, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            iyiz[5] = worky[offsety + 5] * workz[offsetz + 5]; 
            iyiz[6] = worky[offsety + 6] * workz[offsetz + 6]; 
            iyiz[7] = worky[offsety + 7] * workz[offsetz + 7]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data[ijposition] += iyiz[5] * workx[offsetx + 5]; 
                current_data[ijposition] += iyiz[6] * workx[offsetx + 6]; 
                current_data[ijposition] += iyiz[7] * workx[offsetx + 7]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
                current_data_y[ijposition] += iyiz[5] * workx_y[offsetx + 5]; 
                current_data_y[ijposition] += iyiz[6] * workx_y[offsetx + 6]; 
                current_data_y[ijposition] += iyiz[7] * workx_y[offsetx + 7]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR9() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 9 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[9];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 9, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[9];
    double wt[9];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    wt[5] = weights_[offset + 5] * roots_[offset + 5];
    wt[6] = weights_[offset + 6] * roots_[offset + 6];
    wt[7] = weights_[offset + 7] * roots_[offset + 7];
    wt[8] = weights_[offset + 8] * roots_[offset + 8];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    womt[5] = weights_[offset + 5]  - wt[5];
    womt[6] = weights_[offset + 6]  - wt[6];
    womt[7] = weights_[offset + 7]  - wt[7];
    womt[8] = weights_[offset + 8]  - wt[8];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 9, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 9, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            iyiz[5] = worky[offsety + 5] * workz[offsetz + 5]; 
            iyiz[6] = worky[offsety + 6] * workz[offsetz + 6]; 
            iyiz[7] = worky[offsety + 7] * workz[offsetz + 7]; 
            iyiz[8] = worky[offsety + 8] * workz[offsetz + 8]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data[ijposition] += iyiz[5] * workx[offsetx + 5]; 
                current_data[ijposition] += iyiz[6] * workx[offsetx + 6]; 
                current_data[ijposition] += iyiz[7] * workx[offsetx + 7]; 
                current_data[ijposition] += iyiz[8] * workx[offsetx + 8]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
                current_data_y[ijposition] += iyiz[5] * workx_y[offsetx + 5]; 
                current_data_y[ijposition] += iyiz[6] * workx_y[offsetx + 6]; 
                current_data_y[ijposition] += iyiz[7] * workx_y[offsetx + 7]; 
                current_data_y[ijposition] += iyiz[8] * workx_y[offsetx + 8]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR10() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 10 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[10];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 10, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[10];
    double wt[10];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    wt[5] = weights_[offset + 5] * roots_[offset + 5];
    wt[6] = weights_[offset + 6] * roots_[offset + 6];
    wt[7] = weights_[offset + 7] * roots_[offset + 7];
    wt[8] = weights_[offset + 8] * roots_[offset + 8];
    wt[9] = weights_[offset + 9] * roots_[offset + 9];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    womt[5] = weights_[offset + 5]  - wt[5];
    womt[6] = weights_[offset + 6]  - wt[6];
    womt[7] = weights_[offset + 7]  - wt[7];
    womt[8] = weights_[offset + 8]  - wt[8];
    womt[9] = weights_[offset + 9]  - wt[9];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 10, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 10, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            iyiz[5] = worky[offsety + 5] * workz[offsetz + 5]; 
            iyiz[6] = worky[offsety + 6] * workz[offsetz + 6]; 
            iyiz[7] = worky[offsety + 7] * workz[offsetz + 7]; 
            iyiz[8] = worky[offsety + 8] * workz[offsetz + 8]; 
            iyiz[9] = worky[offsety + 9] * workz[offsetz + 9]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data[ijposition] += iyiz[5] * workx[offsetx + 5]; 
                current_data[ijposition] += iyiz[6] * workx[offsetx + 6]; 
                current_data[ijposition] += iyiz[7] * workx[offsetx + 7]; 
                current_data[ijposition] += iyiz[8] * workx[offsetx + 8]; 
                current_data[ijposition] += iyiz[9] * workx[offsetx + 9]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
                current_data_y[ijposition] += iyiz[5] * workx_y[offsetx + 5]; 
                current_data_y[ijposition] += iyiz[6] * workx_y[offsetx + 6]; 
                current_data_y[ijposition] += iyiz[7] * workx_y[offsetx + 7]; 
                current_data_y[ijposition] += iyiz[8] * workx_y[offsetx + 8]; 
                current_data_y[ijposition] += iyiz[9] * workx_y[offsetx + 9]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR11() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 11 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[11];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 11, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[11];
    double wt[11];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    wt[5] = weights_[offset + 5] * roots_[offset + 5];
    wt[6] = weights_[offset + 6] * roots_[offset + 6];
    wt[7] = weights_[offset + 7] * roots_[offset + 7];
    wt[8] = weights_[offset + 8] * roots_[offset + 8];
    wt[9] = weights_[offset + 9] * roots_[offset + 9];
    wt[10] = weights_[offset + 10] * roots_[offset + 10];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    womt[5] = weights_[offset + 5]  - wt[5];
    womt[6] = weights_[offset + 6]  - wt[6];
    womt[7] = weights_[offset + 7]  - wt[7];
    womt[8] = weights_[offset + 8]  - wt[8];
    womt[9] = weights_[offset + 9]  - wt[9];
    womt[10] = weights_[offset + 10]  - wt[10];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 11, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 11, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            iyiz[5] = worky[offsety + 5] * workz[offsetz + 5]; 
            iyiz[6] = worky[offsety + 6] * workz[offsetz + 6]; 
            iyiz[7] = worky[offsety + 7] * workz[offsetz + 7]; 
            iyiz[8] = worky[offsety + 8] * workz[offsetz + 8]; 
            iyiz[9] = worky[offsety + 9] * workz[offsetz + 9]; 
            iyiz[10] = worky[offsety + 10] * workz[offsetz + 10]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data[ijposition] += iyiz[5] * workx[offsetx + 5]; 
                current_data[ijposition] += iyiz[6] * workx[offsetx + 6]; 
                current_data[ijposition] += iyiz[7] * workx[offsetx + 7]; 
                current_data[ijposition] += iyiz[8] * workx[offsetx + 8]; 
                current_data[ijposition] += iyiz[9] * workx[offsetx + 9]; 
                current_data[ijposition] += iyiz[10] * workx[offsetx + 10]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
                current_data_y[ijposition] += iyiz[5] * workx_y[offsetx + 5]; 
                current_data_y[ijposition] += iyiz[6] * workx_y[offsetx + 6]; 
                current_data_y[ijposition] += iyiz[7] * workx_y[offsetx + 7]; 
                current_data_y[ijposition] += iyiz[8] * workx_y[offsetx + 8]; 
                current_data_y[ijposition] += iyiz[9] * workx_y[offsetx + 9]; 
                current_data_y[ijposition] += iyiz[10] * workx_y[offsetx + 10]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR12() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 12 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[12];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 12, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[12];
    double wt[12];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    wt[5] = weights_[offset + 5] * roots_[offset + 5];
    wt[6] = weights_[offset + 6] * roots_[offset + 6];
    wt[7] = weights_[offset + 7] * roots_[offset + 7];
    wt[8] = weights_[offset + 8] * roots_[offset + 8];
    wt[9] = weights_[offset + 9] * roots_[offset + 9];
    wt[10] = weights_[offset + 10] * roots_[offset + 10];
    wt[11] = weights_[offset + 11] * roots_[offset + 11];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    womt[5] = weights_[offset + 5]  - wt[5];
    womt[6] = weights_[offset + 6]  - wt[6];
    womt[7] = weights_[offset + 7]  - wt[7];
    womt[8] = weights_[offset + 8]  - wt[8];
    womt[9] = weights_[offset + 9]  - wt[9];
    womt[10] = weights_[offset + 10]  - wt[10];
    womt[11] = weights_[offset + 11]  - wt[11];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 12, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 12, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            iyiz[5] = worky[offsety + 5] * workz[offsetz + 5]; 
            iyiz[6] = worky[offsety + 6] * workz[offsetz + 6]; 
            iyiz[7] = worky[offsety + 7] * workz[offsetz + 7]; 
            iyiz[8] = worky[offsety + 8] * workz[offsetz + 8]; 
            iyiz[9] = worky[offsety + 9] * workz[offsetz + 9]; 
            iyiz[10] = worky[offsety + 10] * workz[offsetz + 10]; 
            iyiz[11] = worky[offsety + 11] * workz[offsetz + 11]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data[ijposition] += iyiz[5] * workx[offsetx + 5]; 
                current_data[ijposition] += iyiz[6] * workx[offsetx + 6]; 
                current_data[ijposition] += iyiz[7] * workx[offsetx + 7]; 
                current_data[ijposition] += iyiz[8] * workx[offsetx + 8]; 
                current_data[ijposition] += iyiz[9] * workx[offsetx + 9]; 
                current_data[ijposition] += iyiz[10] * workx[offsetx + 10]; 
                current_data[ijposition] += iyiz[11] * workx[offsetx + 11]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
                current_data_y[ijposition] += iyiz[5] * workx_y[offsetx + 5]; 
                current_data_y[ijposition] += iyiz[6] * workx_y[offsetx + 6]; 
                current_data_y[ijposition] += iyiz[7] * workx_y[offsetx + 7]; 
                current_data_y[ijposition] += iyiz[8] * workx_y[offsetx + 8]; 
                current_data_y[ijposition] += iyiz[9] * workx_y[offsetx + 9]; 
                current_data_y[ijposition] += iyiz[10] * workx_y[offsetx + 10]; 
                current_data_y[ijposition] += iyiz[11] * workx_y[offsetx + 11]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}


void SlaterBatch::perform_USVRR13() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 13 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* workx_y = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[13];

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];
    double* current_data_y = &datay_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], 13, worksize, workx, svrr_.svrrfunc[svrr_index]);
    double womt[13];
    double wt[13];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    wt[2] = weights_[offset + 2] * roots_[offset + 2];
    wt[3] = weights_[offset + 3] * roots_[offset + 3];
    wt[4] = weights_[offset + 4] * roots_[offset + 4];
    wt[5] = weights_[offset + 5] * roots_[offset + 5];
    wt[6] = weights_[offset + 6] * roots_[offset + 6];
    wt[7] = weights_[offset + 7] * roots_[offset + 7];
    wt[8] = weights_[offset + 8] * roots_[offset + 8];
    wt[9] = weights_[offset + 9] * roots_[offset + 9];
    wt[10] = weights_[offset + 10] * roots_[offset + 10];
    wt[11] = weights_[offset + 11] * roots_[offset + 11];
    wt[12] = weights_[offset + 12] * roots_[offset + 12];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    womt[2] = weights_[offset + 2]  - wt[2];
    womt[3] = weights_[offset + 3]  - wt[3];
    womt[4] = weights_[offset + 4]  - wt[4];
    womt[5] = weights_[offset + 5]  - wt[5];
    womt[6] = weights_[offset + 6]  - wt[6];
    womt[7] = weights_[offset + 7]  - wt[7];
    womt[8] = weights_[offset + 8]  - wt[8];
    womt[9] = weights_[offset + 9]  - wt[9];
    womt[10] = weights_[offset + 10]  - wt[10];
    womt[11] = weights_[offset + 11]  - wt[11];
    womt[12] = weights_[offset + 12]  - wt[12];
    cix.scale_data_t(workx_y, wt, coeffy_[ii]);
    cix.scale_data(womt, coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], 13, worksize, worky, svrr_.svrrfunc[svrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], 13, worksize, workz, svrr_.svrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
            iyiz[3] = worky[offsety + 3] * workz[offsetz + 3]; 
            iyiz[4] = worky[offsety + 4] * workz[offsetz + 4]; 
            iyiz[5] = worky[offsety + 5] * workz[offsetz + 5]; 
            iyiz[6] = worky[offsety + 6] * workz[offsetz + 6]; 
            iyiz[7] = worky[offsety + 7] * workz[offsetz + 7]; 
            iyiz[8] = worky[offsety + 8] * workz[offsetz + 8]; 
            iyiz[9] = worky[offsety + 9] * workz[offsetz + 9]; 
            iyiz[10] = worky[offsety + 10] * workz[offsetz + 10]; 
            iyiz[11] = worky[offsety + 11] * workz[offsetz + 11]; 
            iyiz[12] = worky[offsety + 12] * workz[offsetz + 12]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                current_data[ijposition] += iyiz[3] * workx[offsetx + 3]; 
                current_data[ijposition] += iyiz[4] * workx[offsetx + 4]; 
                current_data[ijposition] += iyiz[5] * workx[offsetx + 5]; 
                current_data[ijposition] += iyiz[6] * workx[offsetx + 6]; 
                current_data[ijposition] += iyiz[7] * workx[offsetx + 7]; 
                current_data[ijposition] += iyiz[8] * workx[offsetx + 8]; 
                current_data[ijposition] += iyiz[9] * workx[offsetx + 9]; 
                current_data[ijposition] += iyiz[10] * workx[offsetx + 10]; 
                current_data[ijposition] += iyiz[11] * workx[offsetx + 11]; 
                current_data[ijposition] += iyiz[12] * workx[offsetx + 12]; 
                current_data_y[ijposition] = iyiz[0] * workx_y[offsetx + 0]; 
                current_data_y[ijposition] += iyiz[1] * workx_y[offsetx + 1]; 
                current_data_y[ijposition] += iyiz[2] * workx_y[offsetx + 2]; 
                current_data_y[ijposition] += iyiz[3] * workx_y[offsetx + 3]; 
                current_data_y[ijposition] += iyiz[4] * workx_y[offsetx + 4]; 
                current_data_y[ijposition] += iyiz[5] * workx_y[offsetx + 5]; 
                current_data_y[ijposition] += iyiz[6] * workx_y[offsetx + 6]; 
                current_data_y[ijposition] += iyiz[7] * workx_y[offsetx + 7]; 
                current_data_y[ijposition] += iyiz[8] * workx_y[offsetx + 8]; 
                current_data_y[ijposition] += iyiz[9] * workx_y[offsetx + 9]; 
                current_data_y[ijposition] += iyiz[10] * workx_y[offsetx + 10]; 
                current_data_y[ijposition] += iyiz[11] * workx_y[offsetx + 11]; 
                current_data_y[ijposition] += iyiz[12] * workx_y[offsetx + 12]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
  delete[] workx_y;
}



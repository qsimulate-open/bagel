//
// BAGEL - Parallel electron correlation program.
// Filename: bvrr.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/rysint/breitbatch.h>
#include <src/rysint/int2d.h>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace bagel;

#if 0
inline size_t BreitBatch::m(int i, int a, int b, int c, int d) const {
  const int la = basisinfo_[0]->angular_number()+2;
  const int lb = basisinfo_[1]->angular_number()+2;
  const int lc = basisinfo_[2]->angular_number()+2;
  const int ld = basisinfo_[3]->angular_number()+2;
  return i+rank_*(a+la*(b+lb*(c+lc*d)));
}
#endif

void BreitBatch::perform_VRR() {
  const int isize = (amax_+1) * (cmax_+1);
  const int worksize = rank_ * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

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

  const int a = basisinfo_[0]->angular_number();
  const int b = basisinfo_[1]->angular_number();
  const int c = basisinfo_[2]->angular_number();
  const int d = basisinfo_[3]->angular_number();

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;

  const int acsize = size_block_ / primsize_;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * rank_;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D cix(dparamx, &roots_[offset], rank_, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], rank_, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], rank_, worksize, workz, vrr_->vrrfunc[vrr_index]);
  }

  stack_->release(worksize*3, workx);

}

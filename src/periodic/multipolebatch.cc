//
// BAGEL - Parallel electron correlation program.
// Filename: multipolebatch.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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


#include <src/periodic/multipolebatch.h>

using namespace std;
using namespace bagel;

MultipoleBatch::MultipoleBatch(const array<shared_ptr<const Shell>,2>& sh, const shared_ptr<const Atom> atom,
                               shared_ptr<StackMem> stack)
 : MultipoleBatch_base(sh, atom, stack) {

 const double integral_thresh = PRIM_SCREEN_THRESH;
 allocate_arrays(primsize_);

 compute_ss(integral_thresh);
}

void MultipoleBatch::compute() {

  // First, do (a|O_lm|0) using VRR
  const int worksize = num_multipoles_ * (amax_ + 1);

  complex<double>* const workx = stack_->get<complex<double>>(worksize);
  complex<double>* const worky = stack_->get<complex<double>>(worksize);
  complex<double>* const workz = stack_->get<complex<double>>(worksize);

  const int nmul = num_multipoles_;

  for (int index = 0; index != primsize_; ++index) {
    complex<double>* current_data = &data_[index * asize_];
    for (int cnt = 0; cnt != nmul; ++cnt) {
      workx[cnt] = multipole_[index * primsize_ + cnt];
      worky[cnt] = 1.0;
      workz[cnt] = 1.0;
    }

    if (amax_ > 0) {
      int cnt = 0;
      for (int l = 0; l <= ANG_VRR_END; ++l) {
        for (int m = 0; m <= 2*l; ++m, ++cnt) {
          workx[nmul + cnt] = - 2.0 * xb_[index] * AB_[0] * multipole_[index * primsize_ + cnt];
          worky[nmul + cnt] = - 2.0 * xb_[index] * AB_[1] * multipole_[index * primsize_ + cnt];
          workz[nmul + cnt] = - 2.0 * xb_[index] * AB_[2] * multipole_[index * primsize_ + cnt];
          if (abs(m - l) < l - 1) {
            const int i1 = (l - 1) * (l - 1) + m;
            workx[nmul + cnt] += 0.5 * (multipole_[index * primsize_ + i1] + multipole_[index * primsize_ + i1 - 2]);
            worky[nmul + cnt] += complex<double>(0.0, 0.5) * (multipole_[index * primsize_ + i1] + multipole_[index * primsize_ + i1 - 2]);
            workz[nmul + cnt] += multipole_[index * primsize_ + i1 - 1];
          }
        }
      }

      for (int a = 2; a <= amax_; ++a) {
        int cnt = 0;
        for (int l = 0; l <= ANG_VRR_END; ++l) {
          for (int m = 0; m <= 2*l; ++m, ++cnt) {
            workx[a * nmul + cnt] = (a - 1.0) * workx[(a - 2) * nmul + cnt] - 2.0 * xb_[index] * AB_[0] * workx[(a - 1) * nmul + cnt];
            worky[a * nmul + cnt] = (a - 1.0) * worky[(a - 2) * nmul + cnt] - 2.0 * xb_[index] * AB_[1] * worky[(a - 1) * nmul + cnt];
            workz[a * nmul + cnt] = (a - 1.0) * workz[(a - 2) * nmul + cnt] - 2.0 * xb_[index] * AB_[2] * workz[(a - 1) * nmul + cnt];
            if (abs(m - l) < l - 1) {
              const int i1 = (l - 1) * (l - 1) + m;
              workx[a * nmul + cnt] += 0.5 * (workx[(a - 1) * nmul + i1] + multipole_[(a - 1) * nmul + i1 - 2]);
              worky[a * nmul + cnt] += complex<double>(0.0, 0.5) * (worky[(a - 1) * nmul + i1] + multipole_[(a - 1) * nmul + i1 - 2]);
              workz[a * nmul + cnt] += workz[(a - 1) * nmul + i1 - 1];
            }
          }
        }
      }
    }

    // get (a|O|0) from x, y, z components
    const int amax1 = amax_ + 1;
    for (int iz = 0; iz <= amax_; ++iz) {
      for (int iy = 0; iy <= amax_ - iz; ++iy) {
        for (int ix = 0; ix <= amax_ - iy - iz; ++ix) {
          for (int i = 0; i != nmul; ++i) {
            const int pos = ang_mapping_[ix + amax1 * (iy + amax1 * iz)] * nmul + i;
            current_data[pos] = workx[ix * nmul + i] * worky[iy * nmul + i] * workz[iz * nmul + i];
          }
        }
      }
    }
  }

  // contract

  // now get (a|O_lm|b) using HRR


  // release resources
  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void MultipoleBatch::perform_contraction() {

}

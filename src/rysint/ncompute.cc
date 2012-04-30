//
// Newint - Parallel electron correlation program.
// Filename: ncompute.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#define PITWOHALF 17.493418327624862

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <src/rysint/naibatch.h>
#include <src/util/f77.h>
#include <src/rysint/macros.h>
#include <src/rysint/carsphlist.h>
#include <src/stackmem.h>

using namespace std;

extern StackMem* stack;

typedef std::shared_ptr<Atom> RefAtom;

void NAIBatch::compute() {
  const double zero = 0.0;
  const int zeroint = 0;
  const int unit = 1;

  bkup_ = stack->get(size_alloc_);

  const int worksize = rank_ * amax1_;
  
  double* workx = stack->get(worksize);
  double* worky = stack->get(worksize);
  double* workz = stack->get(worksize);

  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);

  double r1x[20];
  double r1y[20];
  double r1z[20];
  double r2[20];

  const int alc = size_alloc_; 
  fill(data_, data_ + alc, zero);

  // perform VRR
  const int natom_unit = natom_ / (2 * L_ + 1);
  assert(natom_ % (2 * L_ + 1) == 0);
  for (int xj = 0; xj != screening_size_; ++xj) {
    const int i = screening_[xj];
    const int iprim = i / natom_;
    const int resid = i % natom_;
    const int cell  = resid / natom_unit - L_;
    const int iatom = resid % natom_unit;
    double disp[3];
    disp[0] = disp[1] = 0.0;
    disp[2] = A_ * cell;
    const int offset_iprim = iprim * asize_;
    double* current_data = &data_[offset_iprim];

    const double* croots = &roots_[i * rank_]; 
    const double* cweights = &weights_[i * rank_]; 
    for (int r = 0; r != rank_; ++r) {
      r1x[r] = p_[i * 3    ] - ax - (p_[i * 3    ] - geom_->atoms(iatom)->position(0) - disp[0]) * croots[r];
      r1y[r] = p_[i * 3 + 1] - ay - (p_[i * 3 + 1] - geom_->atoms(iatom)->position(1) - disp[1]) * croots[r];
      r1z[r] = p_[i * 3 + 2] - az - (p_[i * 3 + 2] - geom_->atoms(iatom)->position(2) - disp[2]) * croots[r];
      r2[r] = (1.0 - croots[r]) * 0.5 / xp_[i];

      workx[r] = cweights[r] * coeff_[i];
      worky[r] = 1.0;
      workz[r] = 1.0;
    }

    if (amax_ > 0) {
      for (int r = 0; r != rank_; ++r) {
        workx[rank_ + r] = r1x[r] * workx[r];
        worky[rank_ + r] = r1y[r] * worky[r];
        workz[rank_ + r] = r1z[r] * workz[r];
      }
      for (int j = 2; j <= amax_; ++j) {
        const int offset = j * rank_;
        for (int r = 0; r != rank_; ++r) {
          workx[offset + r] = r1x[r] * workx[offset - rank_ + r] + (j - 1) * r2[r] * workx[offset - rank_ * 2 + r];
          worky[offset + r] = r1y[r] * worky[offset - rank_ + r] + (j - 1) * r2[r] * worky[offset - rank_ * 2 + r];
          workz[offset + r] = r1z[r] * workz[offset - rank_ + r] + (j - 1) * r2[r] * workz[offset - rank_ * 2 + r];
        }
      }
    }
    // assembly step
    for (int iz = 0; iz <= amax_; ++iz) {
      for (int iy = 0; iy <= amax_ - iz; ++iy) {
        for (int ix = max(0, amin_ - iy - iz); ix <= amax_ - iy - iz; ++ix) {
          const int pos = amapping_[ix + amax1_ * (iy + amax1_ * iz)];
          for (int r = 0; r != rank_; ++r)
          current_data[pos] += workz[iz * rank_ + r] * worky[iy * rank_ + r] * workx[ix * rank_ + r];
        }
      }
    }
  }

  // contract indices 01 
  // data will be stored in bkup_: cont01{ xyz{ } }
  {
    const int m = asize_; 
    perform_contraction(m, data_, prim0size_, prim1size_, bkup_, 
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0size_, 
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1size_);
  }

  // HRR to indices 01
  // data will be stored in data_: cont01{ xyzab{ } }
  {
    if (basisinfo_[1]->angular_number() != 0) { 
      const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      hrr_->hrrfunc_call(hrr_index, contsize_, bkup_, AB_, data_);
    } else {
      ::memcpy(data_, bkup_, size_alloc_ * sizeof(double)); 
    }
  }

  // Cartesian to spherical 01 if necesarry
  // data will be stored in bkup_
  if (spherical_) {
    struct CarSphList carsphlist;
    const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = contsize_;
    carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_); 
  }

  // Sort cont01 and xyzab
  // data will be stored in data_: cont1b{ cont0a{ } }
  if (spherical_) {
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_->sortfunc_call(index, data_, bkup_, cont1size_, cont0size_, 1, swap01_);
  } else {
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_->sortfunc_call(index, bkup_, data_, cont1size_, cont0size_, 1, swap01_);
    ::memcpy(data_, bkup_, size_final_ * sizeof(double)); 
  }

  stack->release(size_alloc_ + rank_ * amax1_ * 3);
}



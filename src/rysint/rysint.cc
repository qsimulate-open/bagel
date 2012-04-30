//
// Newint - Parallel electron correlation program.
// Filename: rysint.cc
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


#include <src/rysint/rysint.h>
#include <cmath>
#include <cassert>

using namespace std;

RysInt::RysInt(const vector<std::shared_ptr<Shell> > info)
 : basisinfo_(info), spherical_(info.front()->spherical()), deriv_rank_(0), tenno_(0),
   hrr_(new HRRList()), vrr_(new VRRList()), sort_(new SortList(info.front()->spherical())) {

}


RysInt::~RysInt() {

}


void RysInt::set_swap_info(const bool swap_bra_ket) {
  // swap 01 indices when needed: Larger angular momentum function comes first
  if (basisinfo_[0]->angular_number() < basisinfo_[1]->angular_number()
   || (basisinfo_[0]->angular_number() == 0 && basisinfo_[1]->angular_number() == 0)) {
    swap(basisinfo_[0], basisinfo_[1]);
    swap01_ = true;
  } else {
    swap01_ = false;
  }
  // swap 23 indices when needed
  if (basisinfo_[2]->angular_number() < basisinfo_[3]->angular_number()
   || (basisinfo_[2]->angular_number() == 0 && basisinfo_[3]->angular_number() == 0)) {
    swap(basisinfo_[2], basisinfo_[3]);
    swap23_ = true;
  } else {
    swap23_ = false;
  }

  if (swap_bra_ket) {
    no_transpose_ = false;
    if (!basisinfo_[0]->angular_number() && !basisinfo_[2]->angular_number()) {
      no_transpose_ = true;
      swap(basisinfo_[0], basisinfo_[2]); 
      swap(basisinfo_[1], basisinfo_[3]); 
      swap(swap01_, swap23_); 
    }
  }
}


void RysInt::set_prim_contsizes() {

}

tuple<int,int,int,int> RysInt::set_angular_info() {
  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  const int ang3 = basisinfo_[3]->angular_number();
  rank_ = ceil(0.5 * (ang0 + ang1 + ang2 + ang3 + 1 + deriv_rank_ + tenno_));
  assert(2 * rank_ >= ang0 + ang1 + ang2 + ang3 + 1 + deriv_rank_ + tenno_); 

  amax_ = ang0 + ang1;
  cmax_ = ang2 + ang3;
  amin_ = ang0;
  cmin_ = ang2;
  amax1_ = amax_ + 1;
  cmax1_ = cmax_ + 1;

  asize_ = 0; 
  for (int i = amin_; i <= amax_; ++i) asize_ += (i + 1) * (i + 2) / 2;
  csize_ = 0; 
  for (int i = cmin_; i <= cmax_; ++i) csize_ += (i + 1) * (i + 2) / 2;

  const int asize_final = (ang0 + 1) * (ang0 + 2) * (ang1 + 1) * (ang1 + 2) / 4;
  const int csize_final = (ang2 + 1) * (ang2 + 2) * (ang3 + 1) * (ang3 + 2) / 4;

  const int asize_final_sph = spherical_ ? (2 * ang0 + 1) * (2 * ang1 + 1) : asize_final;
  const int csize_final_sph = spherical_ ? (2 * ang2 + 1) * (2 * ang3 + 1) : csize_final;

  int cnt = 0;
  for (int i = cmin_; i <= cmax_; ++i) {
    for (int iz = 0; iz <= i; ++iz) { 
      for (int iy = 0; iy <= i - iz; ++iy) { 
        const int ix = i - iy - iz;
        if (ix >= 0) {
          cmapping_[ix + cmax1_ * (iy + cmax1_ * iz)] = cnt;
          ++cnt;
        }
      }
    }
  }
  cnt = 0;
  for (int j = amin_; j <= amax_; ++j) {
    for (int jz = 0; jz <= j; ++jz) { 
      for (int jy = 0; jy <= j - jz; ++jy) { 
        const int jx = j - jy - jz;
        if (jx >= 0){
          amapping_[jx + amax1_ * (jy + amax1_ * jz)] = cnt; 
          ++cnt;
        }
      }
    }
  }
  return make_tuple(asize_final, csize_final, asize_final_sph, csize_final_sph);
}

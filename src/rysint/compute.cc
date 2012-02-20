//
// Newint - Parallel electron correlation program.
// Filename: compute.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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
#include <algorithm>
#include <cstring>
#include <cassert>
#include <src/rysint/eribatch.h>
#include <src/util/f77.h>
#include <src/rysint/macros.h>
#include <src/rysint/carsphlist.h>
#include <src/stackmem.h>

using namespace std;

extern StackMem* stack;

void ERIBatch::compute() {
  bool swapped = false;
  const double zero = 0.0;
  const int zeroint = 0;
  const int unit = 1;

  bkup_ = stack->get(size_alloc_);
  
  const int size = size_alloc_;
  fill(data_, data_ + size, zero);

  // perform VRR
  // data_ will contain the intermediates: prim01{ prim23{ xyz{ } } } 
  switch (rank_) {
    case 1: perform_VRR1(); break;
    case 2: perform_VRR2(); break;
    case 3: perform_VRR3(); break;
    case 4: perform_VRR4(); break;
    case 5: perform_VRR5(); break;
    case 6: perform_VRR6(); break;
    case 7: perform_VRR7(); break;
    case 8: perform_VRR8(); break;
    case 9: perform_VRR9(); break;  
    case 10: perform_VRR10(); break;  
    case 11: perform_VRR11(); break;  
    case 12: perform_VRR12(); break;  
    case 13: perform_VRR13(); break;  
    default: assert(false); break;
  }

  // contract indices 01 
  // data will be stored in bkup_: cont01{ prim23{ xyz{ } } }
  {
    const int m = prim2size_ * prim3size_ * asize_ * csize_; 
    perform_contraction_new_outer(m, data_, prim0size_, prim1size_, bkup_, 
      basisinfo_[0]->contractions(), basisinfo_[0]->contraction_upper(), basisinfo_[0]->contraction_lower(), cont0size_, 
      basisinfo_[1]->contractions(), basisinfo_[1]->contraction_upper(), basisinfo_[1]->contraction_lower(), cont1size_);
  }

  // contract indices 23 
  // data will be stored in data_: cont01{ cont23{ xyz{ } } }
  {
    const int n = cont0size_ * cont1size_;
    perform_contraction_new_inner(n, bkup_, prim2size_, prim3size_, data_, 
      basisinfo_[2]->contractions(), basisinfo_[2]->contraction_upper(), basisinfo_[2]->contraction_lower(), cont2size_, 
      basisinfo_[3]->contractions(), basisinfo_[3]->contraction_upper(), basisinfo_[3]->contraction_lower(), cont3size_);
  }

  // HRR to indices 01
  // data will be stored in bkup_: cont01{ cont23{ xyzf{ xyzab{ } } } }
  if (basisinfo_[1]->angular_number() != 0) { 
    const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    hrr_.hrrfunc_call(hrr_index, contsize_ * csize_, data_, AB_, bkup_);
  } else {
    swapped = (swapped ^ true);
  }

  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  const int ang3 = basisinfo_[3]->angular_number();

  int a = (ang0 + 1) * (ang0 + 2) / 2;
  int b = (ang1 + 1) * (ang1 + 2) / 2;
  int c = (ang2 + 1) * (ang2 + 2) / 2;
  int d = (ang3 + 1) * (ang3 + 2) / 2;
  const int asph = 2 * ang0 + 1;
  const int bsph = 2 * ang1 + 1;
  const int csph = 2 * ang2 + 1;
  const int dsph = 2 * ang3 + 1;


  // Cartesian to spherical 01 if necesarry
  // data will be stored in data_
  struct CarSphList carsphlist;
  const bool need_sph01 = basisinfo_[0]->angular_number() > 1;
  if (spherical_ && need_sph01) {
    const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = contsize_ * csize_;
    if (!swapped)
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_); 
    else
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_); 
    swapped = (swapped ^ true); 
  }


  // Transform batch for the second HRR step
  // data will be stored in data_: cont01{ xyzab{ cont23{ xyzf{ } } } }  if cartesian
  // data will be stored in bkup_: cont01{ xyzab{ cont23{ xyzf{ } } } }  if spherical
  if (basisinfo_[0]->angular_number() != 0) {
    const int m = spherical_ ? (asph * bsph) : (a * b);
    const int n = cont2size_ * cont3size_ * csize_;
    const int nloop = cont0size_ * cont1size_;
    int offset = 0;
    if (swapped) {
      for (int i = 0; i != nloop; ++i, offset += m * n)  
        mytranspose_(&data_[offset], &m, &n, &bkup_[offset]);
    } else {
      for (int i = 0; i != nloop; ++i, offset += m * n)  
        mytranspose_(&bkup_[offset], &m, &n, &data_[offset]);
    }
  } else {
    swapped = (swapped ^ true);
  }

  // HRR to indices 23
  // data will be stored in bkup_: cont01{ xyzab{ cont23{ xyzcd{ } } } } if cartesian
  // data will be stored in data_: cont01{ xyzab{ cont23{ xyzcd{ } } } } if spherical
  if (basisinfo_[3]->angular_number() != 0) { 
    const int hrr_index = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
    if (swapped && spherical_)       hrr_.hrrfunc_call(hrr_index, contsize_ * asph * bsph, bkup_, CD_, data_);
    else if (swapped)                hrr_.hrrfunc_call(hrr_index, contsize_ * a * b, bkup_, CD_, data_);
    else if (!swapped && spherical_) hrr_.hrrfunc_call(hrr_index, contsize_ * asph * bsph, data_, CD_, bkup_);
    else                             hrr_.hrrfunc_call(hrr_index, contsize_ * a * b, data_, CD_, bkup_);
  } else {
    swapped = (swapped ^ true); 
  }

  // Cartesian to spherical 23 if necesarry
  // data will be stored in bkup_
  const bool need_sph23 = basisinfo_[2]->angular_number() > 1;
  if (spherical_ && need_sph23) {
    const int carsphindex = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
    const int nloops = contsize_ * asph * bsph;
    if (swapped)
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_); 
    else
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_); 
    swapped = (swapped ^ true);
  }
  a = asph;
  b = bsph;
  c = csph;
  d = dsph;

  double *data_now = swapped ? bkup_ : data_;
  double *bkup_now = swapped ? data_ : bkup_;

  // Sort cont23 and xyzcd
  // data will be stored in data_: cont01{ xyzab{ cont3d{ cont2c{ } } } }
  if (basisinfo_[2]->angular_number() != 0) {
    const int nloop = a * b * cont0size_ * cont1size_;
    const unsigned int index = basisinfo_[3]->angular_number() * ANG_HRR_END + basisinfo_[2]->angular_number();
    sort_.sortfunc_call(index, data_now, bkup_now, cont3size_, cont2size_, nloop, swap23_);
  } else {
    swapped = (swapped ^ true);
  }

  data_now = swapped ? bkup_ : data_;
  bkup_now = swapped ? data_ : bkup_;
  // transpose batch
  // data will be stored in bkup_: cont3d{ cont2c{ cont01{ xyzab{ } } } } 
  if (!no_transpose_) {
    const int m = c * d * cont2size_ * cont3size_;
    const int n = a * b * cont0size_ * cont1size_; 
    mytranspose_(data_now, &m, &n, bkup_now);
  } else {
    swapped = (swapped ^ true);
  } 

  data_now = swapped ? bkup_ : data_;
  bkup_now = swapped ? data_ : bkup_;
  // Sort cont01 and xyzab
  // data will be stored in data_: cont3d{ cont2c{ cont1b{ cont0a{ } } } }
  if (basisinfo_[0]->angular_number() != 0) {
    const int nloop = c * d * cont2size_ * cont3size_;
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_.sortfunc_call(index, data_now, bkup_now, cont1size_, cont0size_, nloop, swap01_);
  } else {
    swapped = (swapped ^ true);
  }
  
  if (swapped) ::memcpy(data_, bkup_, size * sizeof(double)); 

  stack->release(size_alloc_);
}



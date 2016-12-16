//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: compute.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/integral/carsphlist.h>
#include <src/integral/sortlist.h>
#include <src/integral/hrrlist.h>
#include <src/integral/rys/eribatch.h>

using namespace std;
using namespace bagel;

static const CarSphList carsphlist;
static const HRRList hrr;

void ERIBatch::compute() {
  bool swapped = false;

  double* const stack_save = stack_->get(size_alloc_);
  bkup_ = stack_save;
  fill_n(data_, size_alloc_, 0.0);

  // perform VRR
  // data_ will contain the intermediates: prim01{ prim23{ xyz{ } } }
  switch (rank_) {
    case 1: perform_VRR1(); break;
    case 2: perform_VRR2(); break;
    case 3: perform_VRR3(); break;
    default: perform_VRR(); break;
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
    perform_contraction_new_inner(n, asize_*csize_, bkup_, prim2size_, prim3size_, data_,
      basisinfo_[2]->contractions(), basisinfo_[2]->contraction_upper(), basisinfo_[2]->contraction_lower(), cont2size_,
      basisinfo_[3]->contractions(), basisinfo_[3]->contraction_upper(), basisinfo_[3]->contraction_lower(), cont3size_);
  }

  // HRR to indices 01
  // data will be stored in bkup_: cont01{ cont23{ xyzf{ xyzab{ } } } }
  if (basisinfo_[1]->angular_number() != 0) {
    const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    hrr.hrrfunc_call(hrr_index, contsize_ * csize_, data_, AB_, bkup_);
  } else {
    swapped = swapped != true;
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
  const bool need_sph01 = basisinfo_[0]->angular_number() > 1;
  if (spherical1_ && need_sph01) {
    const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = contsize_ * csize_;
    if (!swapped)
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_);
    else
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_);
    swapped = swapped != true;
    a = asph;
    b = bsph;
  }

  // Transform batch for the second HRR step
  // data will be stored in data_: cont01{ xyzab{ cont23{ xyzf{ } } } }  if cartesian
  // data will be stored in bkup_: cont01{ xyzab{ cont23{ xyzf{ } } } }  if spherical
  if (basisinfo_[0]->angular_number() != 0) {
    const int m = a * b;
    const int n = cont2size_ * cont3size_ * csize_;
    const int nloop = cont0size_ * cont1size_;
    int offset = 0;
    if (swapped) {
      for (int i = 0; i != nloop; ++i, offset += m * n)
        blas::transpose(&data_[offset], m, n, &bkup_[offset]);
    } else {
      for (int i = 0; i != nloop; ++i, offset += m * n)
        blas::transpose(&bkup_[offset], m, n, &data_[offset]);
    }
  } else {
    swapped = swapped != true;
  }

  // HRR to indices 23
  // data will be stored in bkup_: cont01{ xyzab{ cont23{ xyzcd{ } } } } if cartesian
  // data will be stored in data_: cont01{ xyzab{ cont23{ xyzcd{ } } } } if spherical
  if (basisinfo_[3]->angular_number() != 0) {
    const int hrr_index = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
    if (swapped) hrr.hrrfunc_call(hrr_index, contsize_ * a * b, bkup_, CD_, data_);
    else         hrr.hrrfunc_call(hrr_index, contsize_ * a * b, data_, CD_, bkup_);
  } else {
    swapped = swapped != true;
  }

  // Cartesian to spherical 23 if necesarry
  // data will be stored in bkup_
  const bool need_sph23 = basisinfo_[2]->angular_number() > 1;
  if (spherical2_ && need_sph23) {
    const int carsphindex = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
    const int nloops = contsize_ * asph * bsph;
    if (swapped)
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_);
    else
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_);
    swapped = swapped != true;
    c = csph;
    d = dsph;
  }

  // if swapped  bkup contains info
  // if !swapped  data contains info
  double *target_now = swapped ? bkup_ : data_;
  double *source_now = swapped ? data_ : bkup_;

  // Sort cont23 and xyzcd
  // data will be stored in data_: cont01{ xyzab{ cont3d{ cont2c{ } } } }
  if (basisinfo_[2]->angular_number() != 0) {
    const SortList sort2(spherical2_);
    const int nloop = a * b * cont0size_ * cont1size_;
    const unsigned int index = basisinfo_[3]->angular_number() * ANG_HRR_END + basisinfo_[2]->angular_number();
    sort2.sortfunc_call(index, target_now, source_now, cont3size_, cont2size_, nloop, swap23_);
  } else {
    swapped = swapped != true;
  }

  target_now = swapped ? data_ : bkup_;
  source_now = swapped ? bkup_ : data_;
  // transpose batch
  // data will be stored in bkup_: cont3d{ cont2c{ cont01{ xyzab{ } } } }
  if (!swap0123_) {
    const int m = c * d * cont2size_ * cont3size_;
    const int n = a * b * cont0size_ * cont1size_;
    blas::transpose(source_now, m, n, target_now);
  } else {
    swapped = swapped != true;
  }

  target_now = swapped ? bkup_ : data_;
  source_now = swapped ? data_ : bkup_;
  // Sort cont01 and xyzab
  // data will be stored in data_: cont3d{ cont2c{ cont1b{ cont0a{ } } } }
  if (basisinfo_[0]->angular_number() != 0) {
    const SortList sort1(spherical1_);
    const int nloop = c * d * cont2size_ * cont3size_;
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort1.sortfunc_call(index, target_now, source_now, cont1size_, cont0size_, nloop, swap01_);
  } else {
    swapped = swapped != true;
  }

  if (swapped) copy(bkup_, bkup_+size_alloc_, data_);

  stack_->release(size_alloc_, stack_save);
}



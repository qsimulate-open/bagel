//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: scompute.cc
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

#define PITWOHALF 17.493418327624862

#include <src/integral/hrrlist.h>
#include <src/integral/sortlist.h>
#include <src/integral/carsphlist.h>
#include <src/integral/rys/slaterbatch.h>

using namespace std;
using namespace bagel;

#ifdef HAVE_LIBSLATER

const static HRRList hrr;
const static CarSphList carsphlist;

void SlaterBatch::compute() {
  bool swapped = false;

  bkup_ = stack_->get(size_alloc_);
  if (yukawa_)
    bkup2_ = stack_->get(size_alloc_);

  const int size = size_alloc_;
  fill_n(data_, size, 0.0);
  if (yukawa_)
    fill_n(data2_, size, 0.0);

  // perform VRR
  // data_ and data2_ will contain the intermediates: prim01{ prim23{ xyz{ } } }
  if (yukawa_) {
    switch (rank_) {
      case 1: perform_USVRR1(); break;
      case 2: perform_USVRR2(); break;
      default: perform_USVRR(); break;
    }
  } else {
    switch (rank_) {
      case 1: perform_SVRR1(); break;
      case 2: perform_SVRR2(); break;
      default: perform_SVRR(); break;
    }
  }

  // contract indices 01
  // data will be stored in bkup_: cont01{ prim23{ xyz{ } } }
  {
    const int m = prim2size_ * prim3size_ * asize_ * csize_;
    perform_contraction_new_outer(m, data_, prim0size_, prim1size_, bkup_,
            basisinfo_[0]->contractions(), basisinfo_[0]->contraction_upper(), basisinfo_[0]->contraction_lower(), cont0size_,
            basisinfo_[1]->contractions(), basisinfo_[1]->contraction_upper(), basisinfo_[1]->contraction_lower(), cont1size_);
    if (yukawa_)
      perform_contraction_new_outer(m, data2_, prim0size_, prim1size_, bkup2_,
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
    if (yukawa_)
      perform_contraction_new_inner(n, asize_*csize_, bkup2_, prim2size_, prim3size_, data2_,
              basisinfo_[2]->contractions(), basisinfo_[2]->contraction_upper(), basisinfo_[2]->contraction_lower(), cont2size_,
              basisinfo_[3]->contractions(), basisinfo_[3]->contraction_upper(), basisinfo_[3]->contraction_lower(), cont3size_);
  }

  // HRR to indices 01
  // data will be stored in bkup_: cont01{ cont23{ xyzf{ xyzab{ } } } }
  {
    if (basisinfo_[1]->angular_number() != 0) {
      const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      hrr.hrrfunc_call(hrr_index, contsize_ * csize_, data_, AB_, bkup_);
      if (yukawa_)
        hrr.hrrfunc_call(hrr_index, contsize_ * csize_, data2_, AB_, bkup2_);
    } else {
      swapped = true;
    }
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
  if (spherical1_) {
    const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = contsize_ * csize_;
    if (!swapped) {
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_);
      if (yukawa_)
        carsphlist.carsphfunc_call(carsphindex, nloops, bkup2_, data2_);
    } else {
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_);
      if (yukawa_)
        carsphlist.carsphfunc_call(carsphindex, nloops, data2_, bkup2_);
    }
    swapped = (swapped ^ true);
    a = asph;
    b = bsph;
  }


  // Transform batch for the second HRR step
  // data will be stored in data_: cont01{ xyzab{ cont23{ xyzf{ } } } }  if cartesian
  // data will be stored in bkup_: cont01{ xyzab{ cont23{ xyzf{ } } } }  if spherical
  {
    const int m = a * b;
    const int n = cont2size_ * cont3size_ * csize_;
    const int nloop = cont0size_ * cont1size_;
    int offset = 0;
    if (swapped) {
      for (int i = 0; i != nloop; ++i, offset += m * n) {
        blas::transpose(&data_[offset], m, n, &bkup_[offset]);
        if (yukawa_)
          blas::transpose(&data2_[offset], m, n, &bkup2_[offset]);
      }
    } else {
      for (int i = 0; i != nloop; ++i, offset += m * n) {
        blas::transpose(&bkup_[offset], m, n, &data_[offset]);
        if (yukawa_)
          blas::transpose(&bkup2_[offset], m, n, &data2_[offset]);
      }
    }
  }

  // HRR to indices 23
  // data will be stored in bkup_: cont01{ xyzab{ cont23{ xyzcd{ } } } } if cartesian
  // data will be stored in data_: cont01{ xyzab{ cont23{ xyzcd{ } } } } if spherical
  {
    if (basisinfo_[3]->angular_number() != 0) {
      const int hrr_index = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
      if (swapped) hrr.hrrfunc_call(hrr_index, contsize_ * a * b, bkup_, CD_, data_);
      else         hrr.hrrfunc_call(hrr_index, contsize_ * a * b, data_, CD_, bkup_);

      if (yukawa_) {
        if (swapped) hrr.hrrfunc_call(hrr_index, contsize_ * a * b, bkup2_, CD_, data2_);
        else         hrr.hrrfunc_call(hrr_index, contsize_ * a * b, data2_, CD_, bkup2_);
      }
    } else {
      swapped = (swapped ^ true);
    }
  }

  // Cartesian to spherical 23 if necesarry
  // data will be stored in bkup_
  if (spherical2_) {
    const int carsphindex = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
    const int nloops = contsize_ * asph * bsph;
    if (swapped) {
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_);
      if (yukawa_)
        carsphlist.carsphfunc_call(carsphindex, nloops, data2_, bkup2_);
    } else {
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_);
      if (yukawa_)
        carsphlist.carsphfunc_call(carsphindex, nloops, bkup2_, data2_);
    }
    swapped = (swapped ^ true);
    c = csph;
    d = dsph;
  }

  double *data_now = swapped ? bkup_ : data_;
  double *bkup_now = swapped ? data_ : bkup_;
  double *data_now_2 = swapped ? bkup2_ : data2_;
  double *bkup_now_2 = swapped ? data2_ : bkup2_;

  // Sort cont23 and xyzcd
  // data will be stored in data_: cont01{ xyzab{ cont3d{ cont2c{ } } } }
  {
    const SortList sort2(spherical2_);
    const int nloop = a * b * cont0size_ * cont1size_;
    const unsigned int index = basisinfo_[3]->angular_number() * ANG_HRR_END + basisinfo_[2]->angular_number();
    sort2.sortfunc_call(index, data_now, bkup_now, cont3size_, cont2size_, nloop, swap23_);
    if (yukawa_)
      sort2.sortfunc_call(index, data_now_2, bkup_now_2, cont3size_, cont2size_, nloop, swap23_);
  }

  // transpose batch
  // data will be stored in bkup_: cont3d{ cont2c{ cont01{ xyzab{ } } } }
  {
    const int m = c * d * cont2size_ * cont3size_;
    const int n = a * b * cont0size_ * cont1size_;
    blas::transpose(data_now, m, n, bkup_now);
    if (yukawa_)
      blas::transpose(data_now_2, m, n, bkup_now_2);
  }

  // Sort cont01 and xyzab
  // data will be stored in data_: cont3d{ cont2c{ cont1b{ cont0a{ } } } }
  {
    const SortList sort1(spherical1_);
    const int nloop = c * d * cont2size_ * cont3size_;
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort1.sortfunc_call(index, data_now, bkup_now, cont1size_, cont0size_, nloop, swap01_);
    if (yukawa_)
      sort1.sortfunc_call(index, data_now_2, bkup_now_2, cont1size_, cont0size_, nloop, swap01_);
  }

  if (swapped) {
    copy_n(bkup_, size, data_);
    if (yukawa_)
      copy_n(bkup2_, size, data2_);
  }

  if (yukawa_) stack_->release(size_alloc_, bkup2_);
  stack_->release(size_alloc_, bkup_);
}

#endif

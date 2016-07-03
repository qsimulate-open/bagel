//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gcompute.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/util/math/algo.h>
#include <src/integral/sortlist.h>
#include <src/integral/carsphlist.h>
#include <src/integral/rys/gradbatch.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;

void GradBatch::compute() {

  double* const stack_sav = stack_->get(size_block_);
  bkup_ = stack_sav;
  fill_n(data_, size_alloc_, 0.0);
  assert(size_block_*12 == size_alloc_);

  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  const int ang3 = basisinfo_[3]->angular_number();
  const int absize_cart = (ang0+1) * (ang0+2) * (ang1+1) * (ang1+2) / 4;
  const int cdsize_cart = (ang2+1) * (ang2+2) * (ang3+1) * (ang3+2) / 4;

  // perform VRR
  // data_ will contain the intermediates: prim01{ prim23{ xyz{ } } }
  perform_VRR();

  // CAUTION!
  // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order

  // loop over gradient...
  double* cdata = data_;
  for (int i = 0; i != 9; ++i, cdata += size_block_) {
    const int iatom = i/3;
    if (basisinfo_[iatom]->dummy()) continue;

    bool swapped = false;

    double* target_now = bkup_;
    const double* source_now = cdata;
    // contract indices 01
    // data will be stored in bkup_: cont01{ prim23{ xyz{ } } }
    {
      const int m = prim2size_ * prim3size_ * absize_cart * cdsize_cart;
      perform_contraction_new_outer(m, source_now, prim0size_, prim1size_, target_now,
        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_upper(), basisinfo_[0]->contraction_lower(), cont0size_,
        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_upper(), basisinfo_[1]->contraction_lower(), cont1size_);
    }

    target_now = cdata;
    source_now = bkup_;
    // contract indices 23
    // data will be stored in data_: cont01{ cont23{ xyz{ } } }
    {
      const int n = cont0size_ * cont1size_;
      perform_contraction_new_inner(n, absize_cart*cdsize_cart, source_now, prim2size_, prim3size_, target_now,
        basisinfo_[2]->contractions(), basisinfo_[2]->contraction_upper(), basisinfo_[2]->contraction_lower(), cont2size_,
        basisinfo_[3]->contractions(), basisinfo_[3]->contraction_upper(), basisinfo_[3]->contraction_lower(), cont3size_);
    }

    int a = (ang0 + 1) * (ang0 + 2) / 2;
    int b = (ang1 + 1) * (ang1 + 2) / 2;
    int c = (ang2 + 1) * (ang2 + 2) / 2;
    int d = (ang3 + 1) * (ang3 + 2) / 2;
    const int asph = 2 * ang0 + 1;
    const int bsph = 2 * ang1 + 1;
    const int csph = 2 * ang2 + 1;
    const int dsph = 2 * ang3 + 1;

    target_now = bkup_;
    source_now = cdata;
    // Cartesian to spherical 01 if necesarry
    // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order
    const bool need_sph01 = basisinfo_[0]->angular_number() > 1;
    if (spherical1_ && need_sph01) {
      const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      const int nloops = contsize_ * cdsize_cart;
      carsphlist.carsphfunc_call(carsphindex, nloops, source_now, target_now);
      a = asph;
      b = bsph;
    } else {
      swapped = swapped != true;
    }

    target_now = swapped ? bkup_ : cdata;
    source_now = swapped ? cdata : bkup_;
    // the result will be 0(1(x0(x1(2(3(x2(x3)))))))
    if (basisinfo_[0]->angular_number() != 0) {
      const int m = a * b;
      const int n = cont2size_ * cont3size_ * cdsize_cart;
      const int nloop = cont0size_ * cont1size_;
      int offset = 0;
      for (int i = 0; i != nloop; ++i, offset += m * n)
        blas::transpose(source_now+offset, m, n, target_now+offset);
    } else {
      swapped = swapped != true;
    }

    target_now = swapped ? cdata : bkup_;
    source_now = swapped ? bkup_ : cdata;
    // Cartesian to spherical 23 if necesarry
    // data will be stored in bkup_
    // the result will be 0(1(x0(x1(2(3(x2(x3)))))))
    const bool need_sph23 = basisinfo_[2]->angular_number() > 1;
    if (spherical2_ && need_sph23) {
      const int carsphindex = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
      const int nloops = contsize_ * asph * bsph;
      carsphlist.carsphfunc_call(carsphindex, nloops, source_now, target_now);
      c = csph;
      d = dsph;
    } else {
      swapped = swapped != true;
    }

    if (spherical1_) {
    }

    target_now = swapped ? bkup_ : cdata;
    source_now = swapped ? cdata : bkup_;

    // Sort cont23 and xyzcd
    // data will be stored in data_: cont01{ xyzab{ cont3d{ cont2c{ } } } }
    if (basisinfo_[2]->angular_number() != 0) {
      const int nloop = a * b * cont0size_ * cont1size_;
      const unsigned int index = basisinfo_[3]->angular_number() * ANG_HRR_END + basisinfo_[2]->angular_number();
      const SortList sort2(spherical2_);
      sort2.sortfunc_call(index, target_now, source_now, cont3size_, cont2size_, nloop, swap23_);
    } else {
      swapped = swapped != true;
    }

    target_now = swapped ? cdata : bkup_;
    source_now = swapped ? bkup_ : cdata;
    // transpose batch
    // data will be stored in bkup_: cont3d{ cont2c{ cont01{ xyzab{ } } } }
    if (!swap0123_) {
      const int m = c * d * cont2size_ * cont3size_;
      const int n = a * b * cont0size_ * cont1size_;
      blas::transpose(source_now, m, n, target_now);
    } else {
      swapped = swapped != true;
    }

    target_now = swapped ? bkup_ : cdata;
    source_now = swapped ? cdata : bkup_;
    // Sort cont01 and xyzab
    // data will be stored in data_: cont3d{ cont2c{ cont1b{ cont0a{ } } } }
    if (basisinfo_[0]->angular_number() != 0) {
      const int nloop = c * d * cont2size_ * cont3size_;
      const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      const SortList sort1(spherical1_);
      sort1.sortfunc_call(index, target_now, source_now, cont1size_, cont0size_, nloop, swap01_);
    } else {
      swapped = swapped != true;
    }

    if (swapped) copy_n(bkup_, size_block_, cdata);

  } // end of loop 12

  using bagel::blas::ax_plus_y_n;
  if (!basisinfo_[3]->dummy()) {
    // use symmetry for shell3.
    if (!basisinfo_[0]->dummy())
      ax_plus_y_n(-1.0, data_+size_block_*0, 3*size_block_, data_+size_block_*9);
    if (!basisinfo_[1]->dummy())
      ax_plus_y_n(-1.0, data_+size_block_*3, 3*size_block_, data_+size_block_*9);
    if (!basisinfo_[2]->dummy())
      ax_plus_y_n(-1.0, data_+size_block_*6, 3*size_block_, data_+size_block_*9);
  } else {
    // use symmetry for shell2
    assert(!basisinfo_[2]->dummy());
    if (!basisinfo_[0]->dummy())
      ax_plus_y_n(-1.0, data_+size_block_*0, 3*size_block_, data_+size_block_*6);
    if (!basisinfo_[1]->dummy())
      ax_plus_y_n(-1.0, data_+size_block_*3, 3*size_block_, data_+size_block_*6);
  }

  stack_->release(size_block_, stack_sav);
}


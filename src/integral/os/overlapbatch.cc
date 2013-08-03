//
// BAGEL - Parallel electron correlation program.
// Filename: overlapbatch.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#include <cassert>
#include <iostream>
#include <src/integral/hrrlist.h>
#include <src/integral/carsphlist.h>
#include <src/integral/os/overlapbatch.h>
#include <src/util/constants.h>

using namespace std;
using namespace bagel;

const static HRRList hrr;
const static CarSphList carsphlist;


void OverlapBatch::compute() {

  double* const intermediate_p = stack_->get(prim0_ * prim1_ * asize_);
  perform_VRR(intermediate_p);

  double* const intermediate_c = stack_->get(cont0_ * cont1_ * asize_);
  fill(intermediate_c, intermediate_c + cont0_ * cont1_ * asize_, 0.0);
  perform_contraction(asize_, intermediate_p, prim0_, prim1_, intermediate_c,
                      basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                      basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

  double* const intermediate_fi = stack_->get(cont0_ * cont1_ * asize_intermediate_);

  if (basisinfo_[1]->angular_number() != 0) {
    const int hrr_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    hrr.hrrfunc_call(hrr_index, cont0_ * cont1_, intermediate_c, AB_, intermediate_fi);
  } else {
    const unsigned int array_size = cont0_ * cont1_ * asize_intermediate_;
    copy(intermediate_c, intermediate_c+array_size, intermediate_fi);
  }

  if (spherical_) {
    double* const intermediate_i = stack_->get(cont0_ * cont1_ * asize_final_);
    const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = cont0_ * cont1_;
    carsphlist.carsphfunc_call(carsph_index, nloops, intermediate_fi, intermediate_i);

    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_.sortfunc_call(sort_index, data_, intermediate_i, cont1_, cont0_, 1, swap01_);
    stack_->release(cont0_ * cont1_ * asize_final_, intermediate_i);
  } else {
    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_.sortfunc_call(sort_index, data_, intermediate_fi, cont1_, cont0_, 1, swap01_);
  }

  stack_->release(cont0_*cont1_*asize_intermediate_, intermediate_fi);
  stack_->release(cont0_*cont1_*asize_, intermediate_c);
  stack_->release(prim0_*prim1_*asize_, intermediate_p);

}


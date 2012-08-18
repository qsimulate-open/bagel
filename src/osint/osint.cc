//
// Newint - Parallel electron correlation program.
// Filename: osint.cc
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


#include <cmath>
#include <cassert>
#include <iostream>
#include <src/osint/osint.h>
#include <src/grad/goverlapbatch.h>
#include <src/util/constants.h>

using namespace std;


static const double pisqrt__ = ::sqrt(pi__);

OSInt::OSInt(const std::vector<std::shared_ptr<const Shell> >& basis, const int deriv)
 : basisinfo_(basis), spherical_(basis.front()->spherical()), sort_(basis.front()->spherical()), deriv_rank_(deriv), stack_(resources__->get()) {

  assert(basis.size() == 2);

  ang0_ = basisinfo_[0]->angular_number();
  ang1_ = basisinfo_[1]->angular_number();

  if (ang0_ < ang1_) {
    swap(basisinfo_[0], basisinfo_[1]);
    swap(ang0_, ang1_);
    swap01_ = true;
  } else {
    swap01_ = false;
  }

  prim0_ = basisinfo_[0]->num_primitive();
  cont0_ = basisinfo_[0]->num_contracted();
  const vector<double> exponents0 = basisinfo_[0]->exponents();
  prim1_ = basisinfo_[1]->num_primitive();
  cont1_ = basisinfo_[1]->num_contracted();
  const vector<double> exponents1 = basisinfo_[1]->exponents();

  AB_[0] = basisinfo_[0]->position(0) - basisinfo_[1]->position(0);
  AB_[1] = basisinfo_[0]->position(1) - basisinfo_[1]->position(1);
  AB_[2] = basisinfo_[0]->position(2) - basisinfo_[1]->position(2);

  vector<double>::const_iterator expi0, expi1; 
  p_.reserve(3 * prim0_ * prim1_);
  xa_.reserve(prim0_ * prim1_);
  xb_.reserve(prim0_ * prim1_);
  xp_.reserve(prim0_ * prim1_);
  coeffsx_.reserve(prim0_ * prim1_);
  coeffsy_.reserve(prim0_ * prim1_);
  coeffsz_.reserve(prim0_ * prim1_);
  coefftx_.reserve(prim0_ * prim1_);
  coeffty_.reserve(prim0_ * prim1_);
  coefftz_.reserve(prim0_ * prim1_);
  for (expi0 = exponents0.begin(); expi0 != exponents0.end(); ++expi0) {
    for (expi1 = exponents1.begin(); expi1 != exponents1.end(); ++expi1) {
      xa_.push_back(*expi0);
      xb_.push_back(*expi1);
      const double cxp = *expi0 + *expi1;
      const double cxp_inv = 1.0 / cxp; 
      const double px = (basisinfo_[0]->position(0) * *expi0 + basisinfo_[1]->position(0) * *expi1) * cxp_inv;
      const double py = (basisinfo_[0]->position(1) * *expi0 + basisinfo_[1]->position(1) * *expi1) * cxp_inv;
      const double pz = (basisinfo_[0]->position(2) * *expi0 + basisinfo_[1]->position(2) * *expi1) * cxp_inv;
      xp_.push_back(cxp); 
      p_.push_back(px);
      p_.push_back(py);
      p_.push_back(pz);
      const double tmp = pisqrt__ * ::sqrt(cxp_inv); 
      coeffsx_.push_back(tmp * ::exp(- *expi0 * *expi1 * cxp_inv * (AB_[0] * AB_[0])));
      coeffsy_.push_back(tmp * ::exp(- *expi0 * *expi1 * cxp_inv * (AB_[1] * AB_[1])));
      coeffsz_.push_back(tmp * ::exp(- *expi0 * *expi1 * cxp_inv * (AB_[2] * AB_[2])));
      const double xpa = px - basisinfo_[0]->position(0);
      const double ypa = py - basisinfo_[0]->position(1);
      const double zpa = pz - basisinfo_[0]->position(2);
      coefftx_.push_back(coeffsx_.back() * (*expi0 - 2 * *expi0 * *expi0 * (xpa * xpa + 0.5 * cxp_inv)));
      coeffty_.push_back(coeffsy_.back() * (*expi0 - 2 * *expi0 * *expi0 * (ypa * ypa + 0.5 * cxp_inv)));
      coefftz_.push_back(coeffsz_.back() * (*expi0 - 2 * *expi0 * *expi0 * (zpa * zpa + 0.5 * cxp_inv)));
    }
  }

  assert(p_.size() == 3 * prim0_ * prim1_);
  assert(xp_.size() == prim0_ * prim1_);
  assert(coeffsx_.size() == prim0_ * prim1_);
  assert(coefftx_.size() == prim0_ * prim1_);

  amax_ = ang0_ + ang1_ + max(deriv_rank_, 0);
  amax1_ = amax_ + 1;
  amin_ = ang0_;

  asize_ = 0;
  for (int i = amin_; i != amax1_; ++i) asize_ += (i + 1) * (i + 2) / 2;
  asize_intermediate_ = (ang0_ + 1) * (ang0_ + 2) * (ang1_ + 1) * (ang1_ + 2) / 4;
  asize_final_ = spherical_ ? (2 * ang0_ + 1) * (2 * ang1_ + 1) : asize_intermediate_;

  if (deriv_rank_ == 0) {
    size_alloc_ = cont0_ * cont1_ * max(asize_intermediate_, asize_);
    size_block_ = size_alloc_;
  } else if (deriv_rank_ == -1) {
    size_block_ = cont0_ * cont1_ * max(asize_intermediate_, asize_);
    size_alloc_ = size_block_ * 3;
  } else if (deriv_rank_ >= 1) {
    size_alloc_ = prim0_ * prim1_ * asize_intermediate_ * 6; // 3*2
    size_block_ = prim0_ * prim1_ * asize_intermediate_;
  } else {
    throw logic_error("high-order multipoles and derivatives not implemented yet");
  }
  stack_save_ = stack_->get(size_alloc_);
  data_ = stack_save_;

  amapping_.resize(amax1_ * amax1_ * amax1_);
  int cnt = 0;
  for (int i = amin_; i <= amax_; ++i) {
    for (int iz = 0; iz <= i; ++iz) { 
      for (int iy = 0; iy <= i - iz; ++iy) { 
        const int ix = i - iy - iz;
        if (ix >= 0) {
          amapping_[ix + amax1_ * (iy + amax1_ * iz)] = cnt;
          ++cnt;
        }
      }
    }
  }
}

OSInt::~OSInt() {
  stack_->release(size_alloc_, stack_save_);
  resources__->release(stack_);
}


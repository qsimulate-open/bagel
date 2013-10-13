#if 0
//
// BAGEL - Parallel electron correlation program.
// Filename: rysint.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/integral/rys/inline.h>
#include <src/integral/rys/eribatch_base.h>
#include <src/integral/rys/naibatch_base.h>

using namespace bagel;

RysInt::RysInt(const std::array<std::shared_ptr<const Shell>,4>& info, std::shared_ptr<StackMem> stack)
 : basisinfo_(info), spherical1_(info[0]->spherical()), spherical2_(info[2]->spherical()), deriv_rank_(0), tenno_(0), breit_(0), london_(0) {
  assert(spherical1_ == info[1]->spherical());
  assert(spherical2_ == info[3]->spherical());

  if (stack == nullptr) {
    stack_ = resources__->get();
    allocated_here_ = true;
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }
}


RysInt::RysInt(const std::array<std::shared_ptr<const Shell>,2>& info, std::shared_ptr<StackMem> stack)
 : spherical1_(info[0]->spherical()), spherical2_(spherical1_), deriv_rank_(0), tenno_(0), breit_(0) {
  auto dum = std::make_shared<const Shell>(spherical2_);
  basisinfo_ = {{ info[0], info[1], dum, dum }};

  if (stack == nullptr) {
    stack_ = resources__->get();
    allocated_here_ = true;
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }
}


RysInt::~RysInt() {
  // TODO this is a little inconsistent
  // stack should be allocated in the constructor of this class

  stack_->release(size_allocated_, buff_);
  if (tenno_) stack_->release(size_alloc_, stack_save2_);
  stack_->release(size_alloc_, stack_save_);

  if (allocated_here_)
    resources__->release(stack_);
}


void RysInt::set_swap_info(const bool swap_bra_ket) {
  // swap 01 indices when needed: Larger angular momentum function comes first
  if (basisinfo_[0]->angular_number() < basisinfo_[1]->angular_number() || basisinfo_[0]->angular_number() == 0) {
    std::swap(basisinfo_[0], basisinfo_[1]);
    swap01_ = true;
  } else {
    swap01_ = false;
  }
  // swap 23 indices when needed
  if (basisinfo_[2]->angular_number() < basisinfo_[3]->angular_number() || basisinfo_[2]->angular_number() == 0) {
    std::swap(basisinfo_[2], basisinfo_[3]);
    swap23_ = true;
  } else {
    swap23_ = false;
  }

  swap0123_ = false;
  if (swap_bra_ket) {
    if (!basisinfo_[0]->angular_number() && !basisinfo_[2]->angular_number()) {
      swap0123_ = true;
      std::tie(basisinfo_[0], basisinfo_[1], basisinfo_[2], basisinfo_[3], swap01_, swap23_)
        = std::make_tuple(basisinfo_[2], basisinfo_[3], basisinfo_[0], basisinfo_[1], swap23_, swap01_);
      std::swap(spherical1_, spherical2_);
    }
  }
}


void RysInt::set_prim_contsizes() {
  prim0size_ = basisinfo_[0]->num_primitive();
  prim1size_ = basisinfo_[1]->num_primitive();
  prim2size_ = basisinfo_[2]->num_primitive();
  prim3size_ = basisinfo_[3]->num_primitive();
  primsize_ = prim0size_ * prim1size_ * prim2size_ * prim3size_;
  cont0size_ = basisinfo_[0]->num_contracted();
  cont1size_ = basisinfo_[1]->num_contracted();
  cont2size_ = basisinfo_[2]->num_contracted();
  cont3size_ = basisinfo_[3]->num_contracted();
  contsize_ = cont0size_ * cont1size_ * cont2size_ * cont3size_;
}


std::tuple<int,int,int,int> RysInt::set_angular_info() {
  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  const int ang3 = basisinfo_[3]->angular_number();
  rank_ = ceil(0.5 * (ang0 + ang1 + ang2 + ang3 + 1 + deriv_rank_ + tenno_ + breit_));
  assert(2 * rank_ >= ang0 + ang1 + ang2 + ang3 + 1 + deriv_rank_ + tenno_ + breit_);

  amax_ = ang0 + ang1 + deriv_rank_;
  cmax_ = ang2 + ang3 + deriv_rank_;
  amin_ = std::max(ang0 - deriv_rank_, 0);
  cmin_ = std::max(ang2 - deriv_rank_, 0);
  amax1_ = amax_ + 1;
  cmax1_ = cmax_ + 1;

  asize_ = 0;
  csize_ = 0;
  for (int i = amin_; i <= amax_; ++i) asize_ += (i + 1) * (i + 2) / 2;
  for (int i = cmin_; i <= cmax_; ++i) csize_ += (i + 1) * (i + 2) / 2;

  const int asize_final = (ang0 + 1) * (ang0 + 2) * (ang1 + 1) * (ang1 + 2) / 4;
  const int csize_final = (ang2 + 1) * (ang2 + 2) * (ang3 + 1) * (ang3 + 2) / 4;

  const int asize_final_sph = spherical1_ ? (2 * ang0 + 1) * (2 * ang1 + 1) : asize_final;
  const int csize_final_sph = spherical2_ ? (2 * ang2 + 1) * (2 * ang3 + 1) : csize_final;

  int cnt = 0;
  for (int i = cmin_; i <= cmax_; ++i) {
    for (int iz = 0; iz <= i; ++iz) {
      for (int iy = 0; iy <= i - iz; ++iy) {
        const int ix = i - iy - iz;
        if (ix >= 0)
          cmapping_[ix + cmax1_ * (iy + cmax1_ * iz)] = cnt++;
      }
    }
  }
  cnt = 0;
  for (int j = amin_; j <= amax_; ++j) {
    for (int jz = 0; jz <= j; ++jz) {
      for (int jy = 0; jy <= j - jz; ++jy) {
        const int jx = j - jy - jz;
        if (jx >= 0)
          amapping_[jx + amax1_ * (jy + amax1_ * jz)] = cnt++;
      }
    }
  }
  return std::make_tuple(asize_final, csize_final, asize_final_sph, csize_final_sph);
}


void RysInt::set_ab_cd() {
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

  AB_[0] = ax - bx;
  AB_[1] = ay - by;
  AB_[2] = az - bz;
  CD_[0] = cx - dx;
  CD_[1] = cy - dy;
  CD_[2] = cz - dz;
}


void RysInt::allocate_arrays(const size_t ps) {
  size_allocated_ = tenno_ > 0 ? ((rank_ * 2 + 13) * ps) : ((rank_ * 2 + 11) * ps);

  buff_ = stack_->get(size_allocated_);  // stack_->get(size_alloc_) stack_->get((rank_ * 2 + 10) * ps)
  double* pointer = buff_;
  screening_ = (int*)pointer;
                    pointer += ps;
  p_ = pointer;     pointer += ps * 3;
  q_ = pointer;     pointer += ps * 3;
  xp_ = pointer;    pointer += ps;
  xq_ = pointer;    pointer += ps;
  coeff_ = pointer; pointer += ps;
  T_ = pointer;     pointer += ps;
  roots_ = pointer; pointer += rank_ * ps;
  weights_ = pointer; pointer += rank_ * ps;
  if (tenno_) {
    coeffy_ = pointer;pointer += ps;
    U_ = pointer;     pointer += ps;
  }
}


// TODO this is not a good design. Should refactor at certain point using virtual functions...
void RysInt::allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) {
  size_final_ = asize_final_sph * csize_final_sph * contsize_;
  if (deriv_rank_ == 0) {
    const unsigned int size_start = asize_ * csize_ * primsize_;
    const unsigned int size_intermediate = asize_final * csize_ * contsize_;
    const unsigned int size_intermediate2 = asize_final_sph * csize_final * contsize_;
    size_block_ = std::max(size_start, std::max(size_intermediate, size_intermediate2));
    size_alloc_ = size_block_;

    // if this is a two-electron Breit integral
    if (breit_)
      size_alloc_ = 6 * size_block_;

    stack_save_ = stack_->get(size_alloc_);
    stack_save2_ = nullptr;

    // if Slater/Yukawa integrals
    if (tenno_)
      stack_save2_ = stack_->get(size_alloc_);

  // derivative integrals
  } else if (deriv_rank_ == 1) {
    size_block_ = asize_final * csize_final * primsize_;
    // if this is a two-electron gradient integral
    if (dynamic_cast<ERIBatch_base*>(this)) {
      size_alloc_ = 12 * size_block_;
    // if this is an NAI gradient integral
    } else if (dynamic_cast<NAIBatch_base*>(this)) {
      // in this case, we store everything
      size_alloc_ = (dynamic_cast<NAIBatch_base*>(this)->mol()->natom()) * 3.0 * size_block_;
      assert(csize_final == 1);
    } else {
      throw std::logic_error("something is strange in RysInt::allocate_data");
    }
    stack_save_ = stack_->get(size_alloc_);
    stack_save2_ = nullptr;
  }
  data_ = stack_save_;
  data2_ = stack_save2_;
}
#endif

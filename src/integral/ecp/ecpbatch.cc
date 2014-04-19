//
// BAGEL - Parallel electron correlation program.
// Filename: ecpbatch.cc
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

#include <src/integral/ecp/ecpbatch.h>
#include <src/integral/carsphlist.h>
#include <src/integral/sortlist.h>


using namespace bagel;
using namespace std;

const static CarSphList carsphlist;

ECPBatch::ECPBatch(const array<shared_ptr<const Shell>,2>& info, const shared_ptr<const Molecule> mol,
                         shared_ptr<StackMem> stack)
 : basisinfo_(info), mol_(mol) {

  if (stack == nullptr) {
    stack_ = resources__->get();
    allocated_here_ = true;
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }

  integral_thresh_ = PRIM_SCREEN_THRESH;
  max_iter_ = 100;

  spherical_ = basisinfo_[0]->spherical();
  assert(spherical_ == basisinfo_[1]->spherical());

  common_init();

}

void ECPBatch::compute() {

  const double zero = 0.0;
  const SortList sort(spherical_);

  fill_n(data_, size_alloc_, zero);
  double* current_data = &data_[0];
  bkup_ = stack_save_;

  for (auto& aiter : mol_->atoms()) {
    shared_ptr<const ECP> aiter_ecp = aiter->ecp_parameters();

    int index = 0;
    for (int izA = 0; izA <= ang0_; ++izA) {
      for (int iyA = 0; iyA <= ang0_ - izA; ++iyA) {
        const int ixA = ang0_ - izA - iyA;
        const array<int, 3> lA = {izA, iyA, ixA};
        for (int izC = 0; izC <= ang1_; ++izC) {
          for (int iyC = 0; iyC <= ang1_ - izC; ++iyC) {
            const int ixC = ang1_ - izC - iyC;
            const array<int, 3> lC = {izC, iyC, ixC};
            RadialInt<AngularBatch, const shared_ptr<const ECP>, const array<shared_ptr<const Shell>,2>&, const array<int, 3>, const array<int, 3>>
                              radint(aiter_ecp, basisinfo_, lA, lC, false, max_iter_, integral_thresh_);
            current_data[index] = radint.integral();
            ++index;
          }
        }
      }
    }
  }

  {
    perform_contraction(asize_, data_, prim0_, prim1_, bkup_,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);
  }
  copy(bkup_, bkup_ + size_alloc_, data_);

  if (spherical_) {
    const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = cont0_ * cont1_;
    carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_);
  }

  if (spherical_) {
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(index, data_, bkup_, cont1_, cont0_, 1, swap01_);
  } else {
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(index, bkup_, data_, cont1_, cont0_, 1, swap01_);
    copy(bkup_, bkup_ + size_final_, data_);
  }

  stack_->release(size_alloc_, stack_save_);

}

void ECPBatch::perform_contraction(const int asize, const double* prim, const int pdim0, const int pdim1, double* cont,
                           const std::vector<std::vector<double>>& coeff0, const std::vector<std::pair<int, int>>& ranges0, const int cdim0,
                           const std::vector<std::vector<double>>& coeff1, const std::vector<std::pair<int, int>>& ranges1, const int cdim1) {
  const int worksize = pdim1 * asize;
  double* const work = stack_->get<double>(worksize);

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = ranges0[i].first;
    const int end0   = ranges0[i].second;
    std::fill_n(work, worksize, double(0.0));
    for (int j = begin0; j != end0; ++j)
      blas::ax_plus_y_n(coeff0[i][j], prim+j*worksize, worksize, work);

    for (int k = 0; k != cdim1; ++k, cont += asize) {
      const int begin1 = ranges1[k].first;
      const int end1   = ranges1[k].second;
      std::fill_n(cont, asize, double(0.0));
      for (int j = begin1; j != end1; ++j)
        blas::ax_plus_y_n(coeff1[k][j], work+j*asize, asize, cont);
    }
  }
  stack_->release(worksize, work);
}

void ECPBatch::common_init() {

  assert(basisinfo_.size() == 2);

  ang0_ = basisinfo_[0]->angular_number();
  ang1_ = basisinfo_[1]->angular_number();

  if (ang0_ < ang1_) {
    std::swap(basisinfo_[0], basisinfo_[1]);
    std::swap(ang0_, ang1_);
    swap01_ = true;
  } else {
    swap01_ = false;
  }

  prim0_ = basisinfo_[0]->num_primitive();
  cont0_ = basisinfo_[0]->num_contracted();
  prim1_ = basisinfo_[1]->num_primitive();
  cont1_ = basisinfo_[1]->num_contracted();

  amax_ = ang0_ + ang1_;
  amax1_ = amax_ + 1;
  amin_ = ang0_;

  asize_ = (ang0_+1) * (ang0_+2) * (ang1_+1) * (ang1_+2) / 4;
  size_alloc_ = max(asize_ * prim0_ * prim1_, asize_ * cont0_ * cont1_);

  const int asize_sph = spherical_ ? (2 * ang0_ + 1) * (2 * ang1_ + 1) : asize_;
  size_final_ = asize_sph * cont0_ * cont1_;

  stack_save_ = stack_->template get<double>(size_alloc_);
  data_ = stack_save_;

}

ECPBatch::~ECPBatch() {
  if (allocated_here_) resources__->release(stack_);
}

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ecpbatch.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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
  max_iter_ = 20;

  spherical_ = basisinfo_[0]->spherical();
  assert(spherical_ == basisinfo_[1]->spherical());

  common_init();

}


ECPBatch::~ECPBatch() {
  stack_->release(size_alloc_, stack_save_);
  if (allocated_here_) resources__->release(stack_);
}


void ECPBatch::compute() {

  double* const intermediate_c = stack_->get(size_alloc_);
  fill_n(intermediate_c, size_alloc_, 0.0);
  double* const current_data = intermediate_c;

  int i = 0;
  for (int izA = 0; izA <= ang0_; ++izA)
  for (int iyA = 0; iyA <= ang0_ - izA; ++iyA) {
    const int ixA = ang0_ - izA - iyA;
    const array<int, 3> lA = {ixA, iyA, izA};
    for (int izC = 0; izC <= ang1_; ++izC)
    for (int iyC = 0; iyC <= ang1_ - izC; ++iyC) {
      const int ixC = ang1_ - izC - iyC;
      const array<int, 3> lC = {ixC, iyC, izC};
      for (int contA = 0; contA != cont0_; ++contA)
      for (int contC = 0; contC != cont1_; ++contC) {
        double tmp = 0.0;
        for (auto& aiter : mol_->atoms()) {
          shared_ptr<const ECP> aiter_ecp = aiter->ecp_parameters();

          AngularBatch radint(aiter_ecp, basisinfo_, contA, contC, lA, lC, false, max_iter_, integral_thresh_);
          radint.integrate();
          tmp += radint.integral(0);
        }
        const int index = contA * cont1_ * asize_ + contC * asize_ + i;
        current_data[index] = tmp;
      }
      ++i;
    }
  }

  get_data(current_data, data_);

  stack_->release(size_alloc_, intermediate_c);

}

void ECPBatch::get_data(double* intermediate, double* data) {

  fill_n(data, size_alloc_, 0.0);

  double* const intermediate_fi = stack_->get(size_alloc_);
  copy_n(intermediate, size_alloc_, intermediate_fi);

  if (spherical_) {
    double* const intermediate_i = stack_->get(cont0_ * cont1_ * asize_final_);
    const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = cont0_ * cont1_;
    carsphlist.carsphfunc_call(carsph_index, nloops, intermediate_fi, intermediate_i);

    const static SortList sort(true);
    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(sort_index, data, intermediate_i, cont1_, cont0_, 1, swap01_);
    stack_->release(cont0_ * cont1_ * asize_final_, intermediate_i);
  } else {
    const static SortList sort(false);
    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(sort_index, data, intermediate_fi, cont1_, cont0_, 1, swap01_);
  }

  stack_->release(size_alloc_, intermediate_fi);

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

  cont0_ = basisinfo_[0]->num_contracted();
  cont1_ = basisinfo_[1]->num_contracted();

  asize_ = (ang0_+1) * (ang0_+2) * (ang1_+1) * (ang1_+2) / 4;
  asize_final_ = (basisinfo_[0]->spherical() ? (2*ang0_+1) : (ang0_+1)*(ang0_+2)/2)
               * (basisinfo_[1]->spherical() ? (2*ang1_+1) : (ang1_+1)*(ang1_+2)/2);

  size_alloc_ = cont0_ * cont1_ * asize_;

  stack_save_ = stack_->get(size_alloc_);
  data_ = stack_save_;

}

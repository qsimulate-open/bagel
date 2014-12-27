//
// BAGEL - Parallel electron correlation program.
// Filename: multipolebatch_base.cc
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


#include <src/periodic/multipolebatch_base.h>

using namespace std;
using namespace bagel;

MultipoleBatch_base::MultipoleBatch_base(const array<shared_ptr<const Shell>,2>& sh, const shared_ptr<const Atom> a,
                                         shared_ptr<StackMem> stack)
 : basisinfo_(sh), site_(a), spherical_(sh[0]->spherical()) {
  assert(spherical_ = sh[1]->spherical());

  if (stack == nullptr) {
    stack_ = resources__->get();
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }

  init();

}


void MultipoleBatch_base::init() {
  // some definitions for convenience
  AB_[0] = basisinfo_[0]->position(0) - basisinfo_[1]->position(0);
  AB_[1] = basisinfo_[0]->position(1) - basisinfo_[1]->position(1);
  AB_[2] = basisinfo_[0]->position(2) - basisinfo_[1]->position(2);

  // set primitive and contraction sizes
  prim0size_ = basisinfo_[0]->num_primitive();
  prim1size_ = basisinfo_[1]->num_primitive();
  primsize_ = prim0size_ * prim1size_;

  cont0size_ = basisinfo_[0]->num_contracted();
  cont1size_ = basisinfo_[1]->num_contracted();
  contsize_ = cont0size_ * cont1size_;

  // set angular numbers
  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();

  amax_ = ang0 + ang1;
  amax1_ = amax_ + 1;

  asize_ = 0;
  for (int i = 0; i <= amax_; ++i) asize_ += (i + 1) * (i + 2) / 2;

  const int asize_final = (ang0 + 1) * (ang0 + 2) * (ang1 + 1) * (ang1 + 2) / 4;

  const int asize_final_sph = spherical_ ? (2 * ang0 + 1) * (2 * ang1 + 1) : asize_final;

  int cnt = 0;
  for (int i = 0; i <= amax_; ++i) {
    for (int iz = 0; iz <= i; ++iz) {
      for (int iy = 0; iy <= i - iz; ++iy) {
        const int ix = i - iy - iz;
        if (ix >= 0)
          amapping_[ix + amax1_ * (iy + amax1_ * iz)] = cnt++;
      }
    }
  }

  allocate_data(asize_final_sph);
  num_multipoles_ = (ANG_VRR_END + 1) * (ANG_VRR_END + 1);

}


void MultipoleBatch_base::allocate_data(const int asize_final_sph) {

  size_alloc_ = asize_final_sph * contsize_;
}


void MultipoleBatch_base::allocate_arrays(const size_t ps) {

  size_allocated_ = ps * 3;
  buff_ = stack_->get(size_allocated_);
  double* pointer = buff_;
  P_ = pointer;                  pointer += ps * 3;

  multipole_buff = stack_->get<complex<double>>(ps * num_multipoles_);
  complex<double> *cpointer = multipole_buff;
  multipole_ = cpointer;        cpointer += ps * num_multipoles_;
}


void MultipoleBatch_base::compute_ss(const double thr) {

  const vector<double> exp0 = basisinfo_[0]->exponents();
  const vector<double> exp1 = basisinfo_[1]->exponents();

  int index = 0;

  for (auto e0 = exp0.begin(); e0 != exp0.end(); ++e0) {
    for (auto e1 = exp1.begin(); e1 != exp1.end(); ++e1, ++index) {

      const double cxp = *e0 + *e1;
      const double cxp_inv = 1.0 / cxp;
      const double ab = *e0 * *e1;

      P_[index * 3    ] = (basisinfo_[0]->position(0) * *e0 + basisinfo_[1]->position(0) * *e1) * cxp_inv;
      P_[index * 3 + 1] = (basisinfo_[0]->position(1) * *e0 + basisinfo_[1]->position(1) * *e1) * cxp_inv;
      P_[index * 3 + 2] = (basisinfo_[0]->position(2) * *e0 + basisinfo_[1]->position(2) * *e1) * cxp_inv;
      const double Sab = pow(pi__ * cxp_inv, 1.5) * exp(- ab * cxp_inv * (AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2]));

      array<double, 3> PQ;
      PQ[0] = P_[index * 3    ] - site_->position(0);
      PQ[1] = P_[index * 3 + 1] - site_->position(1);
      PQ[2] = P_[index * 3 + 2] - site_->position(2);
      auto Mpq = make_shared<const Multipole>(PQ);
      for (int i = 0; i != num_multipoles_; ++i)
        multipole_[i] = Mpq->multipole(i) * Sab;
    }
  }
}

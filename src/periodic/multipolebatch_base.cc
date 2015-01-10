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
#include <src/periodic/multipole.h>

using namespace std;
using namespace bagel;

MultipoleBatch_base::MultipoleBatch_base(const array<shared_ptr<const Shell>,2>& sh, const array<double, 3> c,
                                         const int lmax, shared_ptr<StackMem> stack)
 : OSIntegral(sh, stack), centre_(c), lmax_(lmax) {

  init();
  common_init();
}


void MultipoleBatch_base::init() {

  assert(max(ang0_, ang1_) < lmax_);
  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1); // spherical
}


void MultipoleBatch_base::allocate_arrays(const size_t ps) {

  multipole_.resize(ps * num_multipoles_);
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

      array<double, 3> P, PQ;
      P[0] = (basisinfo_[0]->position(0) * *e0 + basisinfo_[1]->position(0) * *e1) * cxp_inv;
      P[1] = (basisinfo_[0]->position(1) * *e0 + basisinfo_[1]->position(1) * *e1) * cxp_inv;
      P[2] = (basisinfo_[0]->position(2) * *e0 + basisinfo_[1]->position(2) * *e1) * cxp_inv;
      const double Sab = pow(pi__ * cxp_inv, 1.5) * exp(- ab * cxp_inv * (AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2]));

      PQ[0] = P[0] - centre_[0];
      PQ[1] = P[1] - centre_[1];
      PQ[2] = P[2] - centre_[2];
      auto Mpq = make_shared<const Multipole>(PQ);
      for (int i = 0; i != num_multipoles_; ++i) {
        cout << " size = " << multipole_.size() << " pos = " << index + i * prim0_ * prim1_ << endl;
        multipole_[index + i * prim0_ * prim1_] = Mpq->multipole(i) * Sab;
      }
    }
  }
}

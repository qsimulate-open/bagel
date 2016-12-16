//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexeribatch.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/integral/comprys/complexeribatch.h>
#include <src/integral/comprys/comperirootlist.h>
#include <src/integral/rys/erirootlist.h>
#include <complex>

using namespace std;
using namespace bagel;


ComplexERIBatch::ComplexERIBatch(const array<shared_ptr<const Shell>,4>& _info, const double max_density, const complex<double> dummy, const bool dum,
                   shared_ptr<StackMem> stack)  :  ERIBatch_Base(_info, 0, 0, stack) {

  const double integral_thresh = (max_density != 0.0) ? (PRIM_SCREEN_THRESH / max_density) : 0.0;
  compute_ssss(integral_thresh);

  root_weight(this->primsize_);
}


void ComplexERIBatch::root_weight(const int ps) {
  assert(breit_ == 0);

  // standard 2-D interpolation
  if (any_of(T_, T_+ps, [](complex<double>k){return (k.imag() != 0.0);})
    || any_of(T_, T_+ps, [](complex<double>k){return (k.real() < 0.0);})) {
    complexeriroot__.root(rank_, T_, roots_, weights_, ps);

  // cheaper 1-D interpolation if T_ are all real and positive
  } else {
    const double num = rank_*ps;
    double* const T_temp = stack_->template get<double>(ps);
    double* const roots_temp = stack_->template get<double>(num);
    double* const weights_temp = stack_->template get<double>(num);

    for (int i=0; i!=ps; ++i) {
      T_temp[i] = T_[i].real();
    }
    eriroot__.root(rank_, T_temp, roots_temp, weights_temp, ps);
    for (int i=0; i!=num; ++i) {
      roots_[i] = roots_temp[i];
      weights_[i] = weights_temp[i];
    }

    stack_->release(num, weights_temp);
    stack_->release(num, roots_temp);
    stack_->release(ps, T_temp);
  }
}

complex<double> ComplexERIBatch::get_PQ(const double coord1, const double coord2, const double exp1, const double exp2, const double one12,
                                             const int center1, const int dim, const bool swap) {
  const double Areal = coord1*exp1;
  const double Breal = coord2*exp2;
  const double Aimag = basisinfo_[center1]->vector_potential(dim);
  const double Bimag = basisinfo_[center1+1]->vector_potential(dim);
  const double imag = 0.5 * (swap ? Bimag - Aimag : Aimag - Bimag);
  assert(exp1 != 0.0 || imag == 0.0);
  assert(exp2 != 0.0 || imag == 0.0);
  const complex<double> num (Areal + Breal, imag);
  return num * one12;
}


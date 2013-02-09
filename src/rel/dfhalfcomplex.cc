//
// BAGEL - Parallel electron correlation program.
// Filename: dfhalfcomplex.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/rel/dfhalfcomplex.h>

using namespace std;
using namespace bagel;

DFHalfComplex::DFHalfComplex(shared_ptr<const DFData> df, array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff)
                              : coord_(df->coord()), basis_(df->basis()) {
  alpha_ = df->alpha();
  sigma1_ = shared_ptr<Sigma>(new Sigma(coord_.first));
  sigma2_ = shared_ptr<Sigma>(new Sigma(coord_.second));
  const int index = basis_.first + (coord_.first == Comp::L ? 0 : 2);

  shared_ptr<DFHalfDist> rhalfbj;
  shared_ptr<DFHalfDist> ihalfbj;

  if (df->swapped()) {
    rhalfbj = df->df()->compute_half_transform_swap(rcoeff[index]);
    ihalfbj = df->df()->compute_half_transform_swap(icoeff[index]); 
  } else {
    rhalfbj = df->df()->compute_half_transform(rcoeff[index]);
    ihalfbj = df->df()->compute_half_transform(icoeff[index]); 
  }

  dfhalf_[0] = rhalfbj->apply_J();
  dfhalf_[1] = ihalfbj->apply_J();

  spinor_ = compute_spinor(coord_, basis_);

  shared_ptr<ZMatrix> z1(new ZMatrix(*sigma1_->data()**spinor_.first));
  shared_ptr<ZMatrix> z2(new ZMatrix(*sigma2_->data()**spinor_.second));
  fac_ = (*z1 % *alpha_->data() * *z2).element(0,0);
  assert(fac_ != complex<double>(0.0));

}


void DFHalfComplex::set_sum_diff() {
  df2_[0] = dfhalf_[0]->copy();
  df2_[0]->daxpy(1.0, dfhalf_[1]);
  df2_[1] = dfhalf_[0]->copy();
  df2_[1]->daxpy(-1.0, dfhalf_[1]);
}


void DFHalfComplex::zaxpy(std::complex<double> a, std::shared_ptr<const DFHalfComplex> o) {
  if (imag(a) == 0.0) {
    const double fac = real(a);
    dfhalf_[0]->daxpy(fac, o->dfhalf_[0]);
    dfhalf_[1]->daxpy(fac, o->dfhalf_[1]);
  } else if (real(a) == 0.0) {
    const double fac = imag(a);
    dfhalf_[0]->daxpy(-fac, o->dfhalf_[1]);
    dfhalf_[1]->daxpy( fac, o->dfhalf_[0]);
  } else {
    throw logic_error("DFHalfComplex::zaxpy can be called by real or imaginary coeff (and not complex)");
  }
}

const tuple<int, int> DFHalfComplex::compute_index_Exop(shared_ptr<const DFHalfComplex> o) const {
  const pair<const int, const int> basis2 = o->basis();
  const pair<const int, const int> coord2 = o->coord();

  // 4x4 ZMatrix starting at 0,0 (large, large) or 0,2n (large, small) or 2n,0 (small, large) or 2n,2n (small)
  const int start1 = coord_.first == Comp::L ? 0 : 2;
  const int start2 = coord2.first == Comp::L ? 0 : 2;
  //go from small large to large small or vice versa
  const int index1 = start1 + basis_.second;
  const int index2 = start2 + basis2.second;

  return make_tuple(index1, index2);
}


int DFHalfComplex::coeff_matrix() const {
  return coord_.first == Comp::L ? basis_.second : basis_.second + 2;
}


bool DFHalfComplex::matches(shared_ptr<DFHalfComplex> o) const {
  return coord_.second == o->coord().second && basis_.second == o->basis().second;
}


pair<shared_ptr<ZMatrix>, shared_ptr<ZMatrix>> DFHalfComplex::compute_spinor(pair<const int, const int> coord, pair<const int, const int> basis) {
  pair<shared_ptr<ZMatrix>, shared_ptr<ZMatrix>> spinor;
  spinor.first = shared_ptr<ZMatrix>(new ZMatrix(4,1,true));
  spinor.second = shared_ptr<ZMatrix>(new ZMatrix(4,1,true));
  const int start1 = coord.first == Comp::L ? 0 : 2;
  const int start2 = coord.second == Comp::L ? 0 : 2;
  const int index1 = start1 + basis.first;
  const int index2 = start2 + basis.second;

  spinor.first->element(index1,0) = 1.0;
  spinor.second->element(index2,0) = 1.0;

  return spinor;
}


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

DFHalfComplex::DFHalfComplex(shared_ptr<const DFData> df, std::vector<shared_ptr<ABcases>> bas, array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff)
                              : RelDFBase(*df) {
  common_init();
  basis_ = bas;

  const int index = basis_.front()->basis(0); 
  for (auto& i : basis_)
    if (i->basis(0) != index) throw logic_error("basis should have the same first index");

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

}


DFHalfComplex::DFHalfComplex(array<shared_ptr<DFHalfDist>,2> data, pair<int, int> coord, vector<shared_ptr<ABcases>> bas) : RelDFBase(coord), dfhalf_(data) {
  common_init();
  basis_ = bas; 
}


void DFHalfComplex::set_basis() {
  // does not have to do anything here
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
    const double rfac = real(a);
    dfhalf_[0]->daxpy(rfac, o->dfhalf_[0]);
    dfhalf_[1]->daxpy(rfac, o->dfhalf_[1]);
    const double ifac = imag(a);
    dfhalf_[0]->daxpy(-ifac, o->dfhalf_[1]);
    dfhalf_[1]->daxpy( ifac, o->dfhalf_[0]);
  }
}


bool DFHalfComplex::matches(shared_ptr<DFHalfComplex> o) const {
  return coord_.second == o->coord().second && basis_[0]->basis_second() == o->basis_[0]->basis_second();
}


bool DFHalfComplex::alpha_matches(shared_ptr<DFHalfComplex> o) const {
#if 0
  return basis_[0]->comp() == o->basis_[0]->comp();
#else
  throw logic_error("not yet implemented");
  return true;
#endif
}


bool DFHalfComplex::alpha_matches(shared_ptr<Breit2Index> o) const {
#if 0
  return basis_[0]->comp() == o->comp().second;
#else
  throw logic_error("not yet implemented");
  return true;
#endif
}


shared_ptr<DFHalfComplex> DFHalfComplex::multiply_breit(shared_ptr<Breit2Index> bt) const {
  array<shared_ptr<DFHalfDist>,2> d = {{ dfhalf_[0]->apply_J(bt->k_term()), dfhalf_[1]->apply_J(bt->k_term())}};
  return shared_ptr<DFHalfComplex>(new DFHalfComplex(d, coord_, basis()));
}


list<shared_ptr<DFHalfComplex>> DFHalfComplex::split() {
  list<shared_ptr<DFHalfComplex>> out;
  for (auto i = basis().begin(); i != basis().end(); ++i) {
    if (i == basis().begin()) {
      out.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(dfhalf_, coord_, {*i})));
    } else {
      // TODO Any way to avoid copying?
      array<shared_ptr<DFHalfDist>,2> d = {{ dfhalf_[0]->copy(), dfhalf_[1]->copy() }};
      out.push_back(shared_ptr<DFHalfComplex>(new DFHalfComplex(d, coord_, {*i})));
    }
  }
  return out;
}

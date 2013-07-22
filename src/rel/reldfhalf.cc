//
// BAGEL - Parallel electron correlation program.
// Filename: reldfhalf.cc
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


#include <src/rel/reldfhalf.h>

using namespace std;
using namespace bagel;

RelDFHalf::RelDFHalf(shared_ptr<const RelDF> df, std::vector<shared_ptr<const SpinorInfo>> bas, array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff)
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


RelDFHalf::RelDFHalf(array<shared_ptr<DFHalfDist>,2> data, pair<int, int> coord, vector<shared_ptr<const SpinorInfo>> bas) : RelDFBase(coord), dfhalf_(data) {
  common_init();
  basis_ = bas;
}


RelDFHalf::RelDFHalf(const RelDFHalf& o) : RelDFBase(o.coord_) {
  common_init();
  basis_ = o.basis_;
  dfhalf_[0] = o.dfhalf_[0]->copy();
  dfhalf_[1] = o.dfhalf_[1]->copy();
}


void RelDFHalf::set_sum_diff() {
  df2_[0] = dfhalf_[0]->copy();
  df2_[0]->daxpy(1.0, dfhalf_[1]);
  df2_[1] = dfhalf_[0]->copy();
  df2_[1]->daxpy(-1.0, dfhalf_[1]);
}


void RelDFHalf::zaxpy(std::complex<double> a, std::shared_ptr<const RelDFHalf> o) {
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


bool RelDFHalf::matches(shared_ptr<const RelDFHalf> o) const {
  return coord_.second == o->coord().second && basis_[0]->basis_second() == o->basis_[0]->basis_second() && alpha_matches(o);
}


// WARNING: This Function assumes you have used the split function to make your RelDFHalf object. TODO
bool RelDFHalf::alpha_matches(shared_ptr<const Breit2Index> o) const {
  return basis_[0]->comp() == o->index().second;
}

// WARNING: This Function assumes you have used the split function to make your RelDFHalf objects. TODO
bool RelDFHalf::alpha_matches(shared_ptr<const RelDFHalf> o) const {
  return basis_[0]->comp() == o->basis()[0]->comp();
}


// WARNING: This Function assumes you have used the split function to make your RelDFHalf object. TODO
// Another WARNING: basis(bt) assumes you are multiplying breit2index by 2nd integral, not first
shared_ptr<RelDFHalf> RelDFHalf::multiply_breit2index(shared_ptr<const Breit2Index> bt) const {
  array<shared_ptr<DFHalfDist>,2> d = {{ dfhalf_[0]->apply_J(bt->k_term()), dfhalf_[1]->apply_J(bt->k_term())}};
  return make_shared<RelDFHalf>(d, coord_, new_basis(bt));
}

const vector<shared_ptr<const SpinorInfo>> RelDFHalf::new_basis(shared_ptr<const Breit2Index> bt) const {
  vector<shared_ptr<const SpinorInfo>> out;
  for (auto& i : basis_)
    out.push_back(make_shared<const SpinorInfo>(*i, bt->index().first));
  return out;
}


list<shared_ptr<RelDFHalf>> RelDFHalf::split(const bool docopy) {
  list<shared_ptr<RelDFHalf>> out;
  for (auto i = basis().begin(); i != basis().end(); ++i) {
    if (i == basis().begin() && docopy) {
      out.push_back(make_shared<RelDFHalf>(dfhalf_, coord_, vector<std::shared_ptr<const SpinorInfo>>{*i}));
    } else {
      // TODO Any way to avoid copying?
      array<shared_ptr<DFHalfDist>,2> d = {{ dfhalf_[0]->copy(), dfhalf_[1]->copy() }};
      out.push_back(make_shared<RelDFHalf>(d, coord_, vector<std::shared_ptr<const SpinorInfo>>{*i}));
    }
  }
  return out;
}


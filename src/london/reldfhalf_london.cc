//
// BAGEL - Parallel electron correlation program.
// Filename: reldfhalf_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/london/reldfhalf_london.h>

using namespace std;
using namespace bagel;

RelDFHalf_London::RelDFHalf_London(shared_ptr<const RelDF_London> df, std::vector<shared_ptr<const SpinorInfo>> bas, array<shared_ptr<const Matrix>,4> rcoeff, array<shared_ptr<const Matrix>,4> icoeff)
: RelDFBase(*df) {

  basis_ = bas;

  const int index = basis_.front()->basis(0);
  for (auto& i : basis_)
    if (i->basis(0) != index) throw logic_error("basis should have the same first index");

  // -1 due to dagger (We are transforming the bra.)
  auto icoeff_scaled = make_shared<const Matrix>(*icoeff[index] * (-1.0));

  // TODO Using 4-multiplication - switch to 3
  auto rdf = dynamic_pointer_cast<const DFDist>(df->get_real());
  auto idf = dynamic_pointer_cast<const DFDist>(df->get_imag());
  assert(rdf && idf);
  if (df->swapped()) {
    auto rr = rdf->compute_half_transform_swap(rcoeff[index]);
    auto ri = rdf->compute_half_transform_swap(icoeff_scaled);
    auto ir = idf->compute_half_transform_swap(rcoeff[index]);
    auto ii = idf->compute_half_transform_swap(icoeff_scaled);
    dfhalf_[0] = rr;
    dfhalf_[0]->ax_plus_y(1.0, ii);
    dfhalf_[1] = ri;
    dfhalf_[1]->ax_plus_y(-1.0, ir);
  } else {
    auto rr = rdf->compute_half_transform(rcoeff[index]);
    auto ri = rdf->compute_half_transform(icoeff_scaled);
    auto ir = idf->compute_half_transform(rcoeff[index]);
    auto ii = idf->compute_half_transform(icoeff_scaled);
    dfhalf_[0] = rr;
    dfhalf_[0]->ax_plus_y(-1.0, ii);
    dfhalf_[1] = ri;
    dfhalf_[1]->ax_plus_y(1.0, ir);
  }
}


RelDFHalf_London::RelDFHalf_London(array<shared_ptr<DFHalfDist>,2> data, pair<int, int> cartesian, vector<shared_ptr<const SpinorInfo>> bas) : RelDFBase(cartesian), dfhalf_(data) {
  basis_ = bas;
}


RelDFHalf_London::RelDFHalf_London(const RelDFHalf_London& o) : RelDFBase(o.cartesian_) {
  basis_ = o.basis_;
  dfhalf_[0] = o.dfhalf_[0]->copy();
  dfhalf_[1] = o.dfhalf_[1]->copy();
}


shared_ptr<RelDFHalf_London> RelDFHalf_London::apply_J() const {
  return make_shared<RelDFHalf_London>(array<shared_ptr<DFHalfDist>,2>{{dfhalf_[0]->apply_J(), dfhalf_[1]->apply_J()}}, cartesian_, basis_);
}


void RelDFHalf_London::set_sum_diff() {
  df2_[0] = dfhalf_[0]->copy();
  df2_[0]->ax_plus_y(1.0, dfhalf_[1]);
  df2_[1] = dfhalf_[0]->copy();
  df2_[1]->ax_plus_y(-1.0, dfhalf_[1]);
}


void RelDFHalf_London::ax_plus_y(std::complex<double> a, std::shared_ptr<const RelDFHalf_London> o) {
  if (imag(a) == 0.0) {
    const double fac = real(a);
    dfhalf_[0]->ax_plus_y(fac, o->dfhalf_[0]);
    dfhalf_[1]->ax_plus_y(fac, o->dfhalf_[1]);
  } else if (real(a) == 0.0) {
    const double fac = imag(a);
    dfhalf_[0]->ax_plus_y(-fac, o->dfhalf_[1]);
    dfhalf_[1]->ax_plus_y( fac, o->dfhalf_[0]);
  } else {
    const double rfac = real(a);
    dfhalf_[0]->ax_plus_y(rfac, o->dfhalf_[0]);
    dfhalf_[1]->ax_plus_y(rfac, o->dfhalf_[1]);
    const double ifac = imag(a);
    dfhalf_[0]->ax_plus_y(-ifac, o->dfhalf_[1]);
    dfhalf_[1]->ax_plus_y( ifac, o->dfhalf_[0]);
  }
}


bool RelDFHalf_London::matches(shared_ptr<const RelDFHalf_London> o) const {
  return cartesian_.second == o->cartesian().second && basis_[0]->basis(1) == o->basis_[0]->basis(1) && alpha_matches(o);
}


bool RelDFHalf_London::alpha_matches(shared_ptr<const Breit2Index> o) const {
  assert(basis_.size() == 1);
  return basis_[0]->alpha_comp() == o->index().second;
}


bool RelDFHalf_London::alpha_matches(shared_ptr<const RelDFHalf_London> o) const {
  assert(basis_.size() == 1);
  return basis_[0]->alpha_comp() == o->basis()[0]->alpha_comp();
}


// assumes you are multiplying breit2index by 2nd integral, not first (see alpha_matches)
shared_ptr<RelDFHalf_London> RelDFHalf_London::multiply_breit2index(shared_ptr<const Breit2Index> bt) const {
  throw logic_error("Breit integrals with London orbitals have not yet been implemented.");
  return nullptr;
  /*
  assert(basis_.size() == 1);
  array<shared_ptr<DFHalfDist>,2> d = {{ dfhalf_[0]->apply_J(bt->data()), dfhalf_[1]->apply_J(bt->data())}};

  vector<shared_ptr<const SpinorInfo>> spinor = { make_shared<const SpinorInfo>(basis_[0]->basis(), bt->index().first, bt->index().second) };
  return make_shared<RelDFHalf_London>(d, cartesian_, spinor);
  */
}


list<shared_ptr<RelDFHalf_London>> RelDFHalf_London::split(const bool docopy) {
  list<shared_ptr<RelDFHalf_London>> out;
  for (auto i = basis().begin(); i != basis().end(); ++i) {
    if (i == basis().begin() && docopy) {
      out.push_back(make_shared<RelDFHalf_London>(dfhalf_, cartesian_, vector<std::shared_ptr<const SpinorInfo>>{*i}));
    } else {
      // TODO Any way to avoid copying?
      array<shared_ptr<DFHalfDist>,2> d = {{ dfhalf_[0]->copy(), dfhalf_[1]->copy() }};
      out.push_back(make_shared<RelDFHalf_London>(d, cartesian_, vector<std::shared_ptr<const SpinorInfo>>{*i}));
    }
  }
  return out;
}


shared_ptr<DFDist> RelDFHalfB_London::back_transform(shared_ptr<const Matrix> r, shared_ptr<const Matrix> i, const bool imag) const {
  assert(false);
  shared_ptr<DFDist> out;
  if (!imag) {
    out = dfhalf_[0]->back_transform(r);
    out->ax_plus_y(-1.0, dfhalf_[1]->back_transform(i));
  } else {
    out = dfhalf_[0]->back_transform(i);
    out->ax_plus_y(1.0, dfhalf_[1]->back_transform(r));
  }
  return out;
}

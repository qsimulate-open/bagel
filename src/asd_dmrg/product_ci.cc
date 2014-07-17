//
// BAGEL - Parallel electron correlation program.
// Filename: product_ci.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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


#include <src/asd_dmrg/product_ci.h>

using namespace std;
using namespace bagel;

// Constructor
ProductRASCivec::ProductRASCivec(shared_ptr<RASSpace> space, set<BlockInfo>& left_blocks, const int nelea, const int neleb) :
  space_(space), lblocks_(left_blocks), nelea_(nelea), neleb_(neleb) {

  for (auto& block : lblocks_) {
    const int na = nelea_ - block.nelea;
    const int nb = neleb_ - block.neleb;

    if (na >= 0 && na <= space->norb() && nb >= 0 && nb <= space->norb()) {
      shared_ptr<RASDeterminants> det = space_->det(na, nb);
      sectors_.emplace(block, make_shared<RASBlockVectors>(det, block));
    }
  }
}


/// Copy-constructor
ProductRASCivec::ProductRASCivec(const ProductRASCivec& o) : space_(o.space_), lblocks_(o.lblocks_), nelea_(o.nelea_), neleb_(o.neleb_) {
  for (auto& sec : o.sectors_)
    sectors_.emplace(sec.first, make_shared<RASBlockVectors>(*sec.second));
}


/// Move-constructor
ProductRASCivec::ProductRASCivec(ProductRASCivec&& o) : sectors_(move(o.sectors_)), space_(move(o.space_)),
  lblocks_(move(o.lblocks_)), nelea_(o.nelea_), neleb_(o.neleb_) {}


/// Copy-assignment
ProductRASCivec& ProductRASCivec::operator=(const ProductRASCivec& o) {
  nelea_ = o.nelea_;
  neleb_ = o.neleb_;
  lblocks_ = o.lblocks_;
  space_ = o.space_;
  sectors_.clear();
  for (auto& osec : o.sectors_)
    sectors_.emplace(osec.first, make_shared<RASBlockVectors>(*osec.second));
  return *this;
}


/// Move-assignment
ProductRASCivec& ProductRASCivec::operator=(ProductRASCivec&& o) {
  nelea_ = o.nelea_;
  neleb_ = o.neleb_;
  lblocks_ = move(o.lblocks_);
  space_ = move(o.space_);
  sectors_ = move(o.sectors_);

  return *this;
}


void ProductRASCivec::scale(const double a) {
  for (auto& sec : sectors_) sec.second->scale(a);
}


double ProductRASCivec::dot_product(const ProductRASCivec& o) const {
  assert(matches(o));
  double out = 0.0;
  for (auto& b : lblocks_)
    out += sectors_.at(b)->dot_product(*o.sectors_.at(b));
  return out;
}


void ProductRASCivec::ax_plus_y(const double& a, const ProductRASCivec& o) {
  assert(matches(o));
  for (auto& b : lblocks_)
    sectors_.at(b)->ax_plus_y(a, *o.sectors_.at(b));
}

//
// BAGEL - Parallel electron correlation program.
// Filename: product_civec.cc
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


#include <src/asd_dmrg/product_civec.h>

using namespace std;
using namespace bagel;

// Constructor
ProductRASCivec::ProductRASCivec(shared_ptr<RASSpace> space, shared_ptr<const DMRG_Block> left, const int nelea, const int neleb) :
  space_(space), left_(left), nelea_(nelea), neleb_(neleb) {

  for (auto& block : lblocks()) {
    const int na = nelea_ - block.nelea;
    const int nb = neleb_ - block.neleb;

    if (na >= 0 && na <= space->norb() && nb >= 0 && nb <= space->norb()) {
      shared_ptr<RASDeterminants> det = space_->det(na, nb);
      sectors_.emplace(block, make_shared<RASBlockVectors>(det, block));
    }
  }
}


/// Copy-constructor
ProductRASCivec::ProductRASCivec(const ProductRASCivec& o) : space_(o.space_), left_(o.left_), nelea_(o.nelea_), neleb_(o.neleb_) {
  for (auto& sec : o.sectors_)
    sectors_.emplace(sec.first, make_shared<RASBlockVectors>(*sec.second));
}


/// Move-constructor
ProductRASCivec::ProductRASCivec(ProductRASCivec&& o) : sectors_(move(o.sectors_)), space_(move(o.space_)),
  left_(o.left_), nelea_(o.nelea_), neleb_(o.neleb_) {}


/// Copy-assignment
ProductRASCivec& ProductRASCivec::operator=(const ProductRASCivec& o) {
  nelea_ = o.nelea_;
  neleb_ = o.neleb_;
  left_ = o.left_;
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
  left_ = move(o.left_);
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
  for (auto& b : lblocks())
    out += sectors_.at(b)->dot_product(*o.sectors_.at(b));
  return out;
}


void ProductRASCivec::ax_plus_y(const double& a, const ProductRASCivec& o) {
  assert(matches(o));
  for (auto& b : lblocks())
    sectors_.at(b)->ax_plus_y(a, *o.sectors_.at(b));
}

void ProductRASCivec::print(const double thresh) const {
  for (auto& isec: sectors_) {
    const int nstate = isec.second->nstates();
    for (int ist = 0; ist < nstate; ++ist) {
      cout << "    |na:" << isec.first.nelea << ",nb:" << isec.first.neleb << "," << ist << "> (x)" << endl;
      isec.second->civec(ist).print(thresh);
    }
    cout << endl;
  }
}

// Only the "diagonal" portion is implemented for now
shared_ptr<ProductRASCivec> ProductRASCivec::spin() const {
  auto out = this->clone();
  for (auto& sector : sectors()) {
    // pure block part
    shared_ptr<const Matrix> spinmat = this->left()->spin(sector.first);
    shared_ptr<const RASBlockVectors> source = sector.second;
    shared_ptr<RASBlockVectors> target = out->sector(sector.first);
    dgemm_("N", "T", target->ndim(), target->mdim(), target->mdim(), 1.0, source->data(), source->ndim(), spinmat->data(), spinmat->ndim(),
                                                                     0.0, target->data(), target->ndim());
    // pure ras part
    const int nstate = source->mdim();
    for (int i = 0; i < nstate; ++i) {
      RASCivecView tmpspin = target->civec(i);
      shared_ptr<const RASCivec> sp = source->civec(i).spin();
      tmpspin.ax_plus_y(1.0, *sp);
    }
  }
  return out;
}

double ProductRASCivec::spin_expectation() const { return this->dot_product(*spin()); }

void ProductRASCivec::spin_decontaminate(const double thresh) {
  const int nspin = nelea() - neleb();
  const int max_spin = nelea() + neleb();

  const double pure_expectation = static_cast<double>(nspin * (nspin + 2)) * 0.25;

  shared_ptr<ProductRASCivec> S2 = spin();
  double actual_expectation = dot_product(*S2);

  int k = nspin + 2;
  while( fabs(actual_expectation - pure_expectation) > thresh ) {
    if ( k > max_spin ) { this->print(0.05); throw std::runtime_error("Spin decontamination failed."); }

    const double factor = -4.0/(static_cast<double>(k*(k+2)));
    ax_plus_y(factor, *S2);

    const double norm = this->norm();
    const double rescale = (norm*norm > 1.0e-60) ? 1.0/norm : 0.0;
    scale(rescale);

    S2 = spin();
    actual_expectation = dot_product(*S2);

    k += 2;
  }
}

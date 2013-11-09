//
// BAGEL - Parallel electron correlation program.
// Filename: complexeribatch.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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

#include <src/integral/comprys/complexeribatch.h>
#include <src/integral/comprys/comperirootlist.h>
#include <complex>

using namespace std;
using namespace bagel;

ComplexERIBatch::ComplexERIBatch(const array<shared_ptr<const Shell>,4>& _info, const double max_density, const std::complex<double> dummy, const bool dum,
                   shared_ptr<StackMem> stack) :  ERIBatch_Base(_info, max_density, 0, 0, stack) {

#ifdef LIBINT_INTERFACE
  assert(false);
#endif

  root_weight(this->primsize_);
}

void ComplexERIBatch::root_weight(const int ps) {
  if (breit_ == 0) {
    complexeriroot__.root(rank_, T_, roots_, weights_, ps);
  } else {
    assert(0);
    throw runtime_error("Relativistic calculations have not been set up for London orbitals");
  }
}

void ComplexERIBatch::compute_ssss(const complex<double> integral_thresh) {
  assert(0);
  throw std::runtime_error("compute_ssss not yet implemented for London orbitals");
}


void ComplexERIBatch::compute() {
  assert(0);
  throw logic_error("ComplexERIBatch::compute not yet implemented");
}

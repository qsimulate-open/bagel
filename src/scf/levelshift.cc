//
// BAGEL - Parallel electron correlation program.
// Filename: levelshift.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#include <array>
#include <vector>
#include <string>
#include <memory>
#include <src/scf/levelshift.h>
#include <src/scf/coeff.h>

using namespace bagel;
using namespace std;

ShiftVirtual::ShiftVirtual(const int nocc, const int shift_parameter) : LevelShift(nocc, shift_parameter) {

}

void ShiftVirtual::shift(Matrix1e& coeff) {
  const int nbasis = coeff.geom()->nbasis();

  coeff.add_diag(shift_parameter_, nocc_, nbasis);

  // Done
}

ShiftDimer::ShiftDimer(shared_ptr<const Coeff> monomer_coeff, const double shift_parameter) : LevelShift(monomer_coeff->mdim(), shift_parameter) {
  // Do I actually need to do anything here?
  subspace_ = monomer_coeff;
}

void ShiftDimer::shift(Matrix1e& coeff) {

  
}

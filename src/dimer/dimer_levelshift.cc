//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_levelshift.cc
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

ShiftDimer::ShiftDimer(shared_ptr<const Dimer> dimer, const double shift_parameter) : LevelShift(dimer->nbasis().first, shift_parameter) {
  subspace_ = dimer->proj_coeff()->slice(dimer->nbasis().first, dimer->nbasis().first + dimer->nbasis().second);

  S_ = shared_ptr<const Matrix1e>(new const Overlap(dimer->sgeom()));
  shared_ptr<Matrix1e> S_1(new Matrix1e(S_));
  S_1->inverse();

  subspace_projector_ = shared_ptr<const Matrix1e>(new const Matrix1e( (*subspace_) * (*S_1) * (*(subspace_->transpose())) ));
}

void ShiftDimer::shift(Matrix1e& fock_mo, shared_ptr<const Coeff> coeff) {

  // Project onto subspace
  // Btilde = P_A * B

  // Compute norm
  // norm = Btilde^dagger * S * Btilde

  // Rank according to overlap with subspace
  // Shift the highest
}

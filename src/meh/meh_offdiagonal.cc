//
// BAGEL - Parallel electron correlation program.
// Filename: meh_offdiagonal.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <vector>

#include <src/util/matrix.h>
#include <src/meh/meh.h>

using namespace std;
using namespace bagel;

shared_ptr<Matrix> MultiExcitonHamiltonian::couple_blocks(DimerSubspace& AB, DimerSubspace& ApBp) {
  Coupling term_type = coupling_type(AB, ApBp);

  DimerSubspace* space1 = &AB;
  DimerSubspace* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);

  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    swap(space1,space2);
  }

  shared_ptr<Matrix> out;

  if (term_type == Coupling::none) { out = shared_ptr<Matrix>(new Matrix(space1->dimerstates(), space2->dimerstates())); }
  else if (term_type == Coupling::diagonal) { throw logic_error("Diagonal blocks are being called from the offdiagonal code. Why?"); }
  else if (term_type == Coupling::alphaET) { out = compute_alphaET(*space1, *space2); }
  else if (term_type == Coupling::betaET) { out = compute_betaET(*space1, *space2); }
  else if (term_type == Coupling::ABflip) { out = compute_ABflip(*space1, *space2); }
  else if (term_type == Coupling::alphabetaET) { out = compute_alphabetaET(*space1, *space2); }
  else { throw logic_error("Blocks with an unknown relationship are being coupled"); }

  if (flip) out = out->transpose();

  // fake return statement for now
  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_alphaET(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_betaET(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_ABflip(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_alphabetaET(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  return out;
}

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reloverlap.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <src/util/constants.h>
#include <src/mat1e/rel/reloverlap.h>

using namespace std;
using namespace bagel;

void RelOverlap::compute_() {
  const int n = mol_->nbasis();
  copy_real_block(1.0, 0, 0, n, n, *overlap_);
  copy_real_block(1.0, n, n, n, n, *overlap_);
  copy_real_block(0.5/(c__*c__), 2*n, 2*n, n, n, *kinetic_);
  copy_real_block(0.5/(c__*c__), 3*n, 3*n, n, n, *kinetic_);
}


shared_ptr<ZMatrix> RelOverlap::tildex(const double thresh) const {
  const Matrix tildeo = *overlap_->tildex(thresh);
  Matrix k = tildeo % *kinetic_ * tildeo;
  const bool nosing = k.inverse_half(thresh*1.0e2);
  if (!nosing)
    throw logic_error("positive and negative energy states have different linear dependency");
  const Matrix tildek = tildeo * k;

  const int n = tildeo.ndim();
  const int m = tildeo.mdim();

  auto out = make_shared<ZMatrix>(4*n, 4*m);
  out->copy_real_block(1.0, 0, 0, n, m, tildeo);
  out->copy_real_block(1.0, n, m, n, m, tildeo);
  out->copy_real_block(c__/std::sqrt(0.5), 2*n, 2*m, n, m, tildek);
  out->copy_real_block(c__/std::sqrt(0.5), 3*n, 3*m, n, m, tildek);

  // check numerical stability of the orthogonalization
  assert((*out % *this * *out).is_identity());

  return out;
}


void RelOverlap::inverse() {
#ifndef NDEBUG
  shared_ptr<ZMatrix> ref = this->copy();
#endif
  Matrix oinv(*overlap_);
  oinv.inverse_symmetric();
  Matrix kinv(*kinetic_);
  kinv.inverse_symmetric();
  const int n = oinv.ndim();
  copy_real_block(1.0, 0, 0, n, n, oinv);
  copy_real_block(1.0, n, n, n, n, oinv);
  copy_real_block(2*(c__*c__), 2*n, 2*n, n, n, kinv);
  copy_real_block(2*(c__*c__), 3*n, 3*n, n, n, kinv);

  // check numerical stability of the inversion
#ifndef NDEBUG
  assert((*ref * *this).is_identity());
#endif
}

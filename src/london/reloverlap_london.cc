//
// BAGEL - Parallel electron correlation program.
// Filename: reloverlap_london.cc
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


#include <src/util/constants.h>
#include <src/london/reloverlap_london.h>

using namespace std;
using namespace bagel;

void RelOverlap_London::compute_() {
  const int n = mol_->nbasis();
  auto scalekinetic = make_shared<ZMatrix>(*kinetic_ * (0.5/(c__*c__)));
  copy_block(0, 0, n, n, overlap_);
  copy_block(n, n, n, n, overlap_);
  copy_block(2*n, 2*n, n, n, scalekinetic);
  copy_block(3*n, 3*n, n, n, scalekinetic);
}


shared_ptr<ZMatrix> RelOverlap_London::tildex(const double thresh) const {
  shared_ptr<ZMatrix> tildeo = overlap_->tildex(thresh);
  shared_ptr<ZMatrix> k = make_shared<ZMatrix>(*tildeo % *kinetic_ * *tildeo);
  const bool nosing = k->inverse_half(thresh*1.0e2);
  if (!nosing)
    throw logic_error("positive and negative energy states have different linear dependency");
  shared_ptr<ZMatrix> tildek = make_shared<ZMatrix>(*tildeo * *k);

  *tildek *= (c__/sqrt(0.5));
  const int n = tildeo->ndim();
  const int m = tildeo->mdim();

  auto out = make_shared<ZMatrix>(4*n, 4*m);
  out->copy_block(0, 0, n, m, tildeo);
  out->copy_block(n, m, n, m, tildeo);
  out->copy_block(2*n, 2*m, n, m, tildek);
  out->copy_block(3*n, 3*m, n, m, tildek);
  return out;
}


shared_ptr<ZMatrix> RelOverlap_London::inverse() const {
  shared_ptr<ZMatrix> out = clone();
  ZMatrix oinv(*overlap_);
  oinv.inverse();
  ZMatrix kinv(*kinetic_);
  kinv.inverse();
  kinv *= (2*(c__*c__));
  const int n = oinv.ndim();
  out->copy_block(0, 0, n, n, oinv);
  out->copy_block(n, n, n, n, oinv);
  out->copy_block(2*n, 2*n, n, n, kinv);
  out->copy_block(3*n, 3*n, n, n, kinv);
  return out;
}

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
  ZMatrix scalekinetic = *kinetic_ * (0.5/(c__*c__));
  copy_block(0, 0, n, n, *overlap_);
  copy_block(n, n, n, n, *overlap_);
  copy_block(2*n, 2*n, n, n, scalekinetic);
  copy_block(3*n, 3*n, n, n, scalekinetic);

  const complex<double> r2 (0.25 / (c__*c__));
  const complex<double> i2 (0.0, r2.real());
  add_block( r2*mol_->magnetic_field(2), 2*n, 2*n, n, n, *overlap_);
  add_block(-r2*mol_->magnetic_field(2), 3*n, 3*n, n, n, *overlap_);
  add_block( r2*mol_->magnetic_field(0), 2*n, 3*n, n, n, *overlap_);
  add_block( r2*mol_->magnetic_field(0), 3*n, 2*n, n, n, *overlap_);
  add_block(-i2*mol_->magnetic_field(1), 2*n, 3*n, n, n, *overlap_);
  add_block( i2*mol_->magnetic_field(1), 3*n, 2*n, n, n, *overlap_);
}


shared_ptr<ZMatrix> RelOverlap_London::tildex(const double thresh) const {
  const ZMatrix tildeo = *overlap_->tildex(thresh);

  const int n = tildeo.ndim();
  const int m = tildeo.mdim();

  const int j = mol_->nbasis();
  const ZMatrix soverlap = *get_submatrix(2*j, 2*j, 2*j, 2*j);
  const ZMatrix tildes = *soverlap.tildex(thresh/(c__*c__));
  if (tildes.ndim() != 2*n || tildes.mdim() != 2*m)
    throw logic_error("positive and negative energy states have different linear dependency");

  auto out = make_shared<ZMatrix>(4*n, 4*m);
  out->copy_block(0, 0, n, m, tildeo);
  out->copy_block(n, m, n, m, tildeo);
  out->copy_block(2*n, 2*m, 2*n, 2*m, tildes);
  return out;
}


void RelOverlap_London::inverse() {
  ZMatrix oinv(*overlap_);
  oinv.inverse();
  const int n = oinv.ndim();

  ZMatrix soverlap = *get_submatrix(2*n, 2*n, 2*n, 2*n);
  soverlap.inverse();
  if (soverlap.ndim() != 2*oinv.ndim())
    throw logic_error("positive and negative energy states have different linear dependency");

  copy_block(0, 0, n, n, oinv);
  copy_block(n, n, n, n, oinv);
  copy_block(2*n, 2*n, 2*n, 2*n, soverlap);
}

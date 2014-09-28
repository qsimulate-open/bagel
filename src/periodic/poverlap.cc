//
// BAGEL - Parallel electron correlation program.
// Filename: poverlap.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/poverlap.h>
#include <src/integral/os/overlapbatch.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(POverlap)

POverlap::POverlap(const shared_ptr<const Lattice> lattice) : PMatrix1e(lattice) {

  init(lattice);
}

void POverlap::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Lattice> lattice, const int block) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  OverlapBatch overlap(input);
  overlap.compute();

  data_[block]->copy_block(offsetb1, offsetb0, dimb1, dimb0, overlap.data());

}

shared_ptr<Data> POverlap::tildex(const double thresh_overlap) const {

  auto out = make_shared<Data>(blocksize_, nblock_);
  for (int i = 0; i != nblock_; ++i)
    (*out)[i] = data_[i]->tildex(thresh_overlap);

  return out;
}

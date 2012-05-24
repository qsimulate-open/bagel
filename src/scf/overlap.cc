//
// Newint - Parallel electron correlation program.
// Filename: overlap.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/scf/overlap.h>
#include <src/osint/overlapbatch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>

typedef std::shared_ptr<const Geometry> RefGeometry;
typedef std::shared_ptr<Shell> RefShell;

using namespace std;

Overlap::Overlap(const RefGeometry gm) : Matrix1e(gm) {

  init();
  fill_upper();

}


Overlap::~Overlap() {

}


void Overlap::computebatch(const vector<RefShell>& input, const int offsetb0, const int offsetb1, const int ndim) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis(); 
  const int dimb0 = input[1]->nbasis(); 
  OverlapBatch overlap(input);
  overlap.compute();
  const double* odata = overlap.data();

  int cnt = 0;
  for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
      data_[i * ndim + j] = odata[cnt];
    }
  }
}



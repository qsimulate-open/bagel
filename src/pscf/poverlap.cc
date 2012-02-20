//
// Newint - Parallel electron correlation program.
// Filename: poverlap.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#include <src/pscf/poverlap.h>
#include <src/osint/overlapbatch.h>
#include <iostream>

typedef std::shared_ptr<PGeometry> RefPGeometry;

using namespace std;

POverlap::POverlap(const RefPGeometry g) : PMatrix1e(g) {

  init();

}

POverlap::~POverlap() {

}

void POverlap::computebatch(const vector<shared_ptr<Shell> >& input,
    const int offsetb0, const int offsetb1,
    const int nbasis, const int blockoffset) {

  // input = [b1, b0]
  // because of the convention in integral codes
  assert(input.size() == 2);
  const int dimb0 = input[1]->nbasis(); 
  const int dimb1 = input[0]->nbasis(); 
  OverlapBatch overlap(input);
  overlap.compute();
  const double* odata = overlap.data();

  int cnt = 0;
  for (int j = offsetb0; j != dimb0 + offsetb0; ++j) { // atoms in unit cell 0
    for (int i = offsetb1; i != dimb1 + offsetb1; ++i, ++cnt) { // atoms in unit cell m
      (*data_)[blockoffset + j * nbasis + i] = static_cast<complex<double> >(odata[cnt]);
    }
  }

}


const int POverlap::calculate_thresh() const {

  cout << "  * Maximum overlap" << endl;
  int out = -1;
  for (int i = 0; i <= K(); ++i) {
    double cmax = 0.0;
    const complex<double>* cdata = data_->pointer((i + K()) * blocksize_);
    for (int j = 0; j != blocksize_; ++j, ++cdata) {
      const double current = abs(*cdata);
      cmax = current > cmax ? current : cmax; 
    }
    cout << "   cell " << setw(3) << i << ":  " << cmax << endl;
    if (cmax > 1.0e-15) ++out;
  }
  cout << endl;
  return out;
}

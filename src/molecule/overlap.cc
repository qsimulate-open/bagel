//
// BAGEL - Parallel electron correlation program.
// Filename: overlap.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/molecule/overlap.h>
#include <src/integral/os/overlapbatch.h>

using namespace std;
using namespace bagel;


Overlap::Overlap(const shared_ptr<const Molecule> mo) : Matrix1e(mo) {

  init();
  fill_upper();

}


void Overlap::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();
  OverlapBatch overlap(input);
  overlap.compute();

  copy_block(offsetb1, offsetb0, dimb1, dimb0, overlap.data());
}


shared_ptr<Matrix> Overlap::tildex(const double thresh) const {
  shared_ptr<Matrix> out = this->copy();
  bool nolindep = out->inverse_half(thresh);
  if (!nolindep) {
    // use canonical orthogonalization. Start over
    cout << "    * Using canonical orthogonalization due to linear dependency" << endl << endl;
    out = this->copy();
    unique_ptr<double[]> eig(new double[ndim_]);
    out->diagonalize(eig.get());
    int m = 0;
    for (int i = 0; i != mdim_; ++i) {
      if (eig[i] > thresh) {
        const double e = 1.0/std::sqrt(eig[i]);
        transform(out->element_ptr(0,i), out->element_ptr(0,i+1), out->element_ptr(0,m++), [&e](double a){ return a*e; });
      }
    }
    out = out->slice(0,m);
  }
  return out;
}

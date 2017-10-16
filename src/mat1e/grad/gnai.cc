//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gnai.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#include <src/mat1e/grad/gnai.h>
#include <src/integral/rys/gnaibatch.h>
#include <src/integral/libint/libint.h>

using namespace std;
using namespace bagel;

GNAI::GNAI(shared_ptr<const Molecule> mol) : GMatrix1e(mol) {

  init(mol);
  // fill_upper();

}


void GNAI::computebatch(const array<shared_ptr<const Shell>,2>& s, const vector<int>& a, const vector<int>& o, shared_ptr<const Molecule> m) {

  // s = [b1, b0]
  assert(s.size() == 2);

  GNAIBatch batch(s, m, tie(a[1], a[0]));
  batch.compute();

  const int dimb0 = s[1]->nbasis(), dimb1 = s[0]->nbasis();
  const int natom = m->natom();
  const int offsetb0 = o[0], offsetb1 = o[1];

  for (int l = 0; l != natom; ++l) {
    for (int k = 0; k != 3; ++k) {
      // for (int i = offsetb0, j = 0; i != offsetb0 + dimb0 ; ++i, ++j) {
      //   blas::ax_plus_y_n(1, batch.data(k + 3 * l) + j * dimb1, dimb1, matrices_[k + 3 * l]->element_ptr(offsetb1, i));
      // }
      matrices_[k + 3 * l]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch.data(k + 3 * l));
    }
  }
}



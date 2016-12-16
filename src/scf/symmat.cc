//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: symmat.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/symmat.h>

using namespace std;
using namespace bagel;

SymMat::SymMat(const shared_ptr<const Geometry> gm, const int iop) : Matrix(gm->nbasis(), gm->nbasis()), petite_(gm->plist()) {

  symrot_ = make_shared<SymRotAbel>(petite_->symop(iop), gm->lmax(), gm->spherical());

  const vector<shared_ptr<const Atom>> atoms = gm->atoms();
  const int natom = atoms.size();

  vector<int> offsets(natom, 0);
  for (int i = 1; i != natom; ++i) offsets[i] = offsets[i - 1] + atoms[i - 1]->nbasis();

  const bool spherical = gm->spherical();

  for (int i = 0; i != natom; ++i) {
    const int tatom_num = petite_->sym_atommap(i)[iop];
    const shared_ptr<const Atom> catom = atoms[i];
    const shared_ptr<const Atom> tatom = atoms[tatom_num];

    const vector<shared_ptr<const Shell>> shells = catom->shells();
    const int nshell = shells.size();

    int coffset = offsets[i];
    int toffset = offsets[tatom_num];
    for (int j = 0; j != nshell; ++j) {
      const int ang = shells[j]->angular_number();
      const int size = spherical ? (2 * ang + 1) : ((ang + 1) * (ang + 2) / 2);
      const int nfunc = shells[j]->num_contracted();
      const vector<double> block = symrot_->primrot(ang);
      assert(block.size() == size * size);

      for (int k = 0; k != nfunc; ++k) {
        copy_block(toffset, coffset, size, size, block.data());
        coffset += size;
        toffset += size;
      }
    }
  }

}

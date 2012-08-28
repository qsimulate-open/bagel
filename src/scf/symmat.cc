//
// BAGEL - Parallel electron correlation program.
// Filename: symmat.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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


#include <src/scf/symmat.h>
#include <algorithm>
#include <iostream>
#include <cassert>

using namespace std;
using namespace bagel;

typedef shared_ptr<const Geometry> RefGeometry;
typedef shared_ptr<const Atom> RefAtom;
typedef shared_ptr<const Shell> RefShell;
typedef shared_ptr<Petite> RefPetite;
typedef shared_ptr<SymRotAbel> RefSymRotAbel;

SymMat::SymMat(const RefGeometry gm, const int iop) : Matrix1e(gm), petite_(gm->plist()) {

  symrot_ = RefSymRotAbel(new SymRotAbel(petite_->symop(iop), gm->lmax(), gm->spherical()));

  const vector<RefAtom> atoms = gm->atoms();
  const int natom = atoms.size();
  vector<RefAtom>::const_iterator aiter, biter;

  vector<int> offsets(natom, 0);
  for (int i = 1; i != natom; ++i) offsets[i] = offsets[i - 1] + atoms[i - 1]->nbasis();

  const bool spherical = gm->spherical();

  for (int i = 0; i != natom; ++i) {
    const int tatom_num = petite_->sym_atommap(i)[iop];
    const RefAtom catom = atoms[i];
    const RefAtom tatom = atoms[tatom_num];

    const vector<RefShell> shells = catom->shells();
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
        int index = 0;
        for (int ic = coffset; ic != coffset + size; ++ic) {
          for (int it = toffset; it != toffset + size; ++it, ++index)
            data_[ic * nbasis_ + it] = block[index];
        }
        coffset += size;
        toffset += size;
      }
    }
  }

}


SymMat::~SymMat() {

}

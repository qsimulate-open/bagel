//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: small1e.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/mat1e/rel/small1e.h>

using namespace std;
using namespace bagel;

template<> void Small1e<NAIBatch>::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule> mol) {
  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  if (mol->natom() < nucleus_blocksize__) {
    SmallInts1e<NAIBatch, shared_ptr<const Molecule>> batch(input, mol);
    batch.compute();
    for (int i = 0; i != this->Nblocks(); ++i)
      this->matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch[i]);
  } else {
    const vector<shared_ptr<const Molecule>> atom_subsets = mol->split_atoms(nucleus_blocksize__);
    for (auto& current_mol : atom_subsets) {
      SmallInts1e<NAIBatch, shared_ptr<const Molecule>> batch(input, current_mol);
      batch.compute();
      for (int i = 0; i != this->Nblocks(); ++i)
        this->matrices_[i]->add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, batch[i]);
    }
  }

}


template<> void Small1e<ERIBatch>::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule> mol) {
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();
#ifdef LIBINT_INTERFACE
  SmallInts1e<Libint, shared_ptr<const Molecule>> batch(input, mol);
#else
  SmallInts1e<ERIBatch, shared_ptr<const Molecule>> batch(input, mol);
#endif

  vector<shared_ptr<const Shell>> nshells;
  for (auto& i : mol->atoms()) {
    if (i->finite_nucleus()) {
      const double fac = - i->atom_charge()*pow(i->atom_exponent()/pi__, 1.5);
      nshells.push_back(make_shared<Shell>(i->spherical(), i->position(), 0, vector<double>{i->atom_exponent()},
                        vector<vector<double>>{{fac}}, vector<pair<int,int>>{make_pair(0,1)}));
    }
  }
  batch.compute(nshells);

  for (int i = 0; i != this->Nblocks(); ++i)
    this->matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch[i]);
}


//
// BAGEL - Parallel electron correlation program.
// Filename: small1e_london.cc
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

#include <src/london/small1e_london.h>

using namespace std;
using namespace bagel;

template<> void Small1e_London<ComplexERIBatch>::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, std::shared_ptr<const Molecule> mol) {
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();
  SmallInts1e_London<ComplexERIBatch> batch(input, mol);

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


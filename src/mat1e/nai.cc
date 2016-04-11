//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: nai.cc
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


#include <src/mat1e/nai.h>
#include <src/integral/rys/naibatch.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/libint/libint.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(NAI)

NAI::NAI(shared_ptr<const Molecule> mol) : Matrix1e(mol) {

  init(mol);
  fill_upper();

}


void NAI::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule> mol) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  {
    NAIBatch nai(input, mol);
    nai.compute();
    copy_block(offsetb1, offsetb0, dimb1, dimb0, nai.data());
  }

  if (mol->has_finite_nucleus()) {
    auto dummy = make_shared<const Shell>(input[0]->spherical());
    for (auto& i : mol->atoms()) {
      if (i->finite_nucleus()) {
        const double fac = - i->atom_charge()*pow(i->atom_exponent()/pi__, 1.5);
        auto in = make_shared<Shell>(i->spherical(), i->position(), 0, vector<double>{i->atom_exponent()}, vector<vector<double>>{{fac}}, vector<pair<int,int>>{make_pair(0,1)});
        const array<shared_ptr<const Shell>,4> shells{{ dummy, in, input[0], input[1] }};
#ifdef LIBINT_INTERFACE
        Libint eri(shells);
#else
        ERIBatch eri(shells, 2.0);
#endif
        eri.compute();
        add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, eri.data());
      }
    }
  }

}



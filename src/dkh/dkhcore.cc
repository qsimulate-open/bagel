//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcore.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Raymond Wang <yiqunwang2021@u.northwestern.edu> 
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


#include <src/dkh/dkhcore.h>
#include <src/integral/os/kineticbatch.h>
#include <src/integral/os/mmbatch.h>
#include <src/integral/rys/naibatch.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/rys/r0batch.h>
#include <src/integral/rys/r1batch.h>
#include <src/integral/rys/r2batch.h>
#include <src/integral/libint/libint.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(DKHcore)

DKHcore::DKHcore(shared_ptr<const Molecule> mol) : Matrix1e(mol) {

  init(mol);
  cout << "Using DKH" << endl;
  fill_upper();
}

void DKHcore::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule> mol) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  {
    KineticBatch kinetic(input);
    kinetic.compute();

    copy_block(offsetb1, offsetb0, dimb1, dimb0, kinetic.data());
  }
  
  {
    if (mol->natom() < nucleus_blocksize__) {
      NAIBatch nai(input, mol);
      nai.compute();
      add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, nai.data());
    } else {
      const vector<shared_ptr<const Molecule>> atom_subsets = mol->split_atoms(nucleus_blocksize__);
      for (auto& current_mol : atom_subsets) {
        NAIBatch nai(input, current_mol);
        nai.compute();
        add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, nai.data());
      }
    }
  }

}



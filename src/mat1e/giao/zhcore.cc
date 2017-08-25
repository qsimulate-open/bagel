//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zhcore.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/mat1e/giao/zhcore.h>
#include <src/integral/compos/complexkineticbatch.h>
#include <src/integral/comprys/complexnaibatch.h>
#include <src/integral/comprys/complexeribatch.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(ZHcore)

ZHcore::ZHcore(shared_ptr<const Molecule> mol, shared_ptr<const HcoreInfo> hcoreinfo/*nothing here for now*/) : ZMatrix1e(mol) {

  init(mol);
  fill_upper_conjg();

  if (mol->atoms().front()->use_ecp_basis())
    throw runtime_error("ECP is not available with a GIAO basis.");
}


void ZHcore::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule> mol) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  {
    ComplexKineticBatch kinetic(input, mol->magnetic_field());
    kinetic.compute();

    copy_block(offsetb1, offsetb0, dimb1, dimb0, kinetic.data());
  }
  {
    if (mol->natom() < nucleus_blocksize__) {
      ComplexNAIBatch nai(input, mol);
      nai.compute();
      add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, nai.data());
    } else {
      const vector<shared_ptr<const Molecule>> atom_subsets = mol->split_atoms(nucleus_blocksize__);
      for (auto& current_mol : atom_subsets) {
        ComplexNAIBatch nai(input, current_mol);
        nai.compute();
        add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, nai.data());
      }
    }
  }

  if (mol->has_finite_nucleus()) {
    auto dummy = make_shared<const Shell>(input[0]->spherical());
    for (auto& i : mol->atoms()) {
      if (i->finite_nucleus()) {
        const double fac = - i->atom_charge()*pow(i->atom_exponent()/pi__, 1.5);
        auto in = make_shared<Shell>(i->spherical(), i->position(), 0, vector<double>{i->atom_exponent()}, vector<vector<double>>{{fac}}, vector<pair<int,int>>{make_pair(0,1)});
        const array<shared_ptr<const Shell>,4> shells{{ dummy, in, input[0], input[1] }};
        ComplexERIBatch eri(shells, 2.0);
        eri.compute();
        add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, eri.data());
      }
    }
  }

  if (mol->external()) throw std::logic_error("Wasn't planning to compute dipole in an external electric field with London orbitals");
#if 0
    DipoleBatch dipole(input, mol);
    dipole.compute();
    const size_t block = dipole.size_block();
    const double* dip = dipole.data();

    int cnt = 0;
    for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
      for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
        data_[i*ndim_ + j] += dip[cnt        ]*mol->external(0);
        data_[i*ndim_ + j] += dip[cnt+block  ]*mol->external(1);
        data_[i*ndim_ + j] += dip[cnt+block*2]*mol->external(2);
      }
    }
#endif
}



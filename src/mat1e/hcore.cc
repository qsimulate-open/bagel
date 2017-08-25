//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hcore.cc
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


#include <src/mat1e/hcore.h>
#include <src/integral/os/kineticbatch.h>
#include <src/integral/os/mmbatch.h>
#include <src/integral/rys/naibatch.h>
#include <src/integral/rys/eribatch.h>
#include <src/integral/rys/r0batch.h>
#include <src/integral/rys/r1batch.h>
#include <src/integral/rys/r2batch.h>
#include <src/integral/ecp/ecpbatch.h>
#include <src/integral/ecp/soecpbatch.h>
#include <src/integral/libint/libint.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(Hcore)

Hcore::Hcore(shared_ptr<const Molecule> mol, shared_ptr<const HcoreInfo> hcoreinfo) : Matrix1e(mol), hso_(make_shared<HSO>(mol->nbasis())) {
  if (hcoreinfo->standard() || hcoreinfo->ecp()) {
    init(mol);
    fill_upper();
  } else {
    auto hcore = hcoreinfo->compute(mol);
    copy_n(hcore->data(), hcore->size(), data());
  }
}

void Hcore::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule> mol) {

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

  if (mol->atoms().front()->use_ecp_basis()) {
    {
      R0Batch r0(input, mol);
      r0.compute();

      add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, r0.data());
    }
    {
      R1Batch r1(input, mol);
      r1.compute();

      add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, r1.data());
    }
    {
      R2Batch r2(input, mol);
      r2.compute();

      add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, r2.data());
    }
    {
      ECPBatch ecp(input, mol);
      ecp.compute();

      add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, ecp.data());
    }
    {
      SOECPBatch soecp(input, mol);
      soecp.compute();

      hso_->construct_iaa(offsetb1, offsetb0, dimb1, dimb0, soecp.data());
      hso_->construct_rab(offsetb1, offsetb0, dimb1, dimb0, soecp.data1());
      hso_->construct_iab(offsetb1, offsetb0, dimb1, dimb0, soecp.data2());
    }
  }

  if (mol->has_finite_nucleus()) {
    auto dummy = make_shared<const Shell>(input[0]->spherical());
    for (auto& i : mol->atoms()) {
      if (i->finite_nucleus()) {
        const double fac = - i->atom_charge()*pow(i->atom_exponent()/pi__, 1.5);
        auto in = make_shared<Shell>(i->spherical(), i->position(), 0, vector<double>{i->atom_exponent()}, vector<vector<double>>{{fac}}, vector<pair<int,int>>{{0,1}});
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

  if (mol->external()) {
    DipoleBatch dipole(input, mol);
    dipole.compute();
    const size_t block = dipole.size_block();
    const double* dip = dipole.data();

    int cnt = 0;
    for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
      for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
        element(j, i) += dip[cnt        ]*mol->external(0);
        element(j, i) += dip[cnt+block  ]*mol->external(1);
        element(j, i) += dip[cnt+block*2]*mol->external(2);
      }
    }
  }
}



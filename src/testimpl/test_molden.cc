//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_molden.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker < shane.parker@u.northwestern.edu >
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


#include <src/util/io/moldenout.h>
#include <src/util/io/moldenin.h>
#include <src/scf/hf/fock.h>

double molden_out_energy(std::string inp1, std::string inp2) {

  auto ofs = std::make_shared<std::ofstream>(inp1 + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  {
    std::stringstream ss; ss << location__ << inp1 << ".json";
    auto idata = std::make_shared<const PTree>(ss.str());
    auto keys = idata->get_child("bagel");
    std::shared_ptr<Geometry> geom;

    std::shared_ptr<const Reference> ref;

    for (auto& itree : *keys) {
      const std::string method = to_lower(itree->get<std::string>("title", ""));

      if (method == "molecule") {
        geom = std::make_shared<Geometry>(itree);

      } else if (method == "hf") {
        auto scf = std::make_shared<RHF>(itree, geom);
        scf->compute();
        ref = scf->conv_to_ref();
      }
      else if (method == "print") {
        std::string out_file = itree->get<std::string>("file", inp1 + ".molden");

        MoldenOut mfs(out_file);
        mfs << geom;
        mfs << ref;

      }
    }
  }

  double energy = 0.0;
  {
    std::stringstream ss; ss << location__ << inp2 << ".json";
    auto idata = std::make_shared<const PTree>(ss.str());
    std::shared_ptr<const PTree> keys = idata->get_child("bagel");
    std::shared_ptr<const PTree> mol = *keys->begin();

    const std::string method = to_lower(mol->get<std::string>("title", ""));
    if (method != "molecule") throw std::logic_error("broken test case");
    auto geom = std::make_shared<Geometry>(mol);

    std::string filename = inp1 + ".molden";
    MoldenIn mfs(filename, geom->spherical());
    mfs.read();
    MOInfo mo(geom);
    mfs >> mo;
    std::shared_ptr<const Coeff> coeff = mo.coeff;

    std::shared_ptr<const Matrix> ao_density = coeff->form_density_rhf(geom->nele()/2);
    auto hcore = std::make_shared<const Hcore>(geom);
    auto fock = std::make_shared<const Fock<1>>(geom, hcore, ao_density, geom->schwarz());

    auto hcore_fock = std::make_shared<const Matrix>(*hcore + *fock);
    energy = ((*ao_density)*(*hcore_fock)).trace();
    energy = 0.5*energy + geom->nuclear_repulsion();

  }
  std::cout.rdbuf(backup_stream);
  return energy;
}

BOOST_AUTO_TEST_SUITE(TEST_MOLDEN)

BOOST_AUTO_TEST_CASE(MOLDEN) {
    BOOST_CHECK(compare(molden_out_energy("hf_write_mol_sph", "hf_read_mol_sph"),   -99.84772354 ));
    BOOST_CHECK(compare(molden_out_energy("hf_write_mol_cart", "hf_read_mol_cart"), -99.84911270 ));
}

BOOST_AUTO_TEST_SUITE_END()

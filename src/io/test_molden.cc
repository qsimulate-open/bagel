//
// BAGEL - Parallel electron correlation program.
// Filename: test_molden.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker < shane.parker@u.northwestern.edu >
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


#include <src/io/moldenout.h>
#include <src/io/moldenin.h>
#include <src/scf/fock.h>

double molden_out_energy(std::string inp1, std::string inp2) {

    std::shared_ptr<std::ofstream> ofs(new std::ofstream(inp1 + ".testout", std::ios::trunc));
    std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  {
    std::shared_ptr<InputData> idata(new InputData("../../test/" + inp1 + ".in"));
    std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
    std::list<std::pair<std::string, std::multimap<std::string, std::string>>> keys = idata->data();

    std::shared_ptr<Reference> ref;

    for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
      if (iter->first == "df-hf") {
        std::shared_ptr<SCF<1>> scf(new SCF<1>(iter->second, geom));
        scf->compute();
        ref = scf->conv_to_ref();
      }
      else if (iter->first == "print") {
        std::multimap<std::string, std::string> pdata = iter->second;
        bool orbitals = read_input<bool>(pdata, "orbitals", false);
        std::string out_file = read_input<std::string>(pdata, "file", inp1 + ".molden");
     
        MoldenOut mfs(out_file);
        mfs << geom;
        mfs << ref;
        mfs.close();

      }
    }
  }

  double energy = 0.0;
  {
    std::shared_ptr<InputData> idata(new InputData("../../test/" + inp2 + ".in"));
    std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
    std::list<std::pair<std::string, std::multimap<std::string, std::string>>> keys = idata->data();

    std::shared_ptr<const Coeff> coeff(new Coeff(geom));

    std::string filename = inp1 + ".molden";
    MoldenIn mfs(filename, geom->spherical());
    mfs.read();
    mfs >> coeff;

    std::shared_ptr<Matrix> ao_density = coeff->form_density_rhf(geom->nele()/2);
    std::shared_ptr<const Matrix> hcore(new Hcore(geom));
    std::shared_ptr<Fock<1>> fock(new Fock<1>(geom, hcore, ao_density, geom->schwarz()));

    std::shared_ptr<Matrix> hcore_fock(new Matrix(*hcore + *fock));
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

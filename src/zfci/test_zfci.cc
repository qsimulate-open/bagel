//
// BAGEL - Parallel electron correlation program.
// Filename: test_zfci.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/zfci/zharrison.h>

std::vector<double> relfci_energy(std::string inp) {

  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = "../../test/" + inp + ".json";
  auto idata = std::make_shared<const PTree>(filename);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<const Reference> ref;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "dhf") {
      auto scf = std::make_shared<Dirac>(itree, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "zfci") {
      auto fci = std::make_shared<ZHarrison>(itree, geom, ref);
      fci->compute();
      std::cout.rdbuf(backup_stream);
      return fci->energy();
    }
  }
  assert(false);
  return std::vector<double>();
}

std::vector<double> reference_relfci_energy() {
  std::vector<double> out(1);
  out[0] = -98.63281881;
  return out;
}

std::vector<double> reference_relfci_energy2() {
  std::vector<double> out(1);
  out[0] = -98.62253505;
  return out;
}

std::vector<double> reference_relfci_energy3() {
  std::vector<double> out(1);
  out[0] = -98.62292308;
  return out;
}


BOOST_AUTO_TEST_SUITE(TEST_RELFCI)

BOOST_AUTO_TEST_CASE(ZHARRISON) {
  BOOST_CHECK(compare(relfci_energy("hf_sto3g_relfci_coulomb"), reference_relfci_energy()));
  BOOST_CHECK(compare(relfci_energy("hf_sto3g_relfci_gaunt"), reference_relfci_energy2()));
  BOOST_CHECK(compare(relfci_energy("hf_sto3g_relfci_breit"), reference_relfci_energy3()));
}

BOOST_AUTO_TEST_SUITE_END()

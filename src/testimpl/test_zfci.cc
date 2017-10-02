//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_zfci.cc
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


#include <src/ci/zfci/relfci.h>

std::vector<double> relfci_energy(std::string inp) {

  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = location__ + inp + ".json";
  auto idata = std::make_shared<const PTree>(filename);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<const Reference> ref;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "dhf") {
      auto scf = std::make_shared<Dirac>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "zfci") {
      auto fci = std::make_shared<RelFCI>(itree, geom, ref);
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
  out[0] = -98.62253450;
  return out;
}

std::vector<double> reference_relfci_energy3() {
  std::vector<double> out(1);
  out[0] = -98.62292281;
  return out;
}

std::vector<double> reference_relfci_energy4() {
  std::vector<double> out(2);
  out[0] = -3189.56011628;
  out[1] = -3185.82830884;
  return out;
}

std::vector<double> reference_relfci_energy5() {
  std::vector<double> out(1);
  out[0] = -1.28898724;
  return out;
}

std::vector<double> reference_relfci_energy6() {
  std::vector<double> out(2);
  out[0] = -679.32861057;
  out[1] = -679.23201607;
  return out;
}

std::vector<double> reference_relfci_energy7() {
  std::vector<double> out(1);
  out[0] = -100.06081004;
  return out;
}

BOOST_AUTO_TEST_SUITE(TEST_RELFCI)

BOOST_AUTO_TEST_CASE(ZHARRISON) {
  BOOST_CHECK(compare(relfci_energy("hf_sto3g_relfci_coulomb"), reference_relfci_energy()));
  BOOST_CHECK(compare(relfci_energy("hf_sto3g_relfci_gaunt"), reference_relfci_energy2()));
  BOOST_CHECK(compare(relfci_energy("hf_sto3g_relfci_breit"), reference_relfci_energy3()));
}

BOOST_AUTO_TEST_SUITE_END()

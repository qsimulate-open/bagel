//
// BAGEL - Parallel electron correlation program.
// Filename: test_ras.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <src/ras/rasci.h>

std::vector<double> ras_energy(std::string inp) {

  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = "../../test/" + inp + ".bgl";
  auto idata = std::make_shared<const PTree>(filename);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<const Reference> ref;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<SCF>(itree, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "rohf"){
      auto scf = std::make_shared<ROHF>(itree, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "ras") {
      std::shared_ptr<RASCI> ras = std::make_shared<RASCI>(itree, geom, ref);

      ras->compute();
      std::cout.rdbuf(backup_stream);
      return ras->energy();
    }
  }
  assert(false);
  return std::vector<double>();
}

std::vector<double> reference_ras_energy_h2o_full() {
  std::vector<double> out(2);
  out[0] = -75.28693755;
  out[1] = -74.86783699;
  return out;
}
std::vector<double> reference_ras_energy_hhe_full() {
  std::vector<double> out(2);
  out[0] = -3.30315793;
  out[1] = -2.78407879;
  return out;
}
std::vector<double> reference_ras_energy_h2o_restricted() {
  std::vector<double> out(2);
  out[0] = -75.28621267;
  out[1] = -74.84580572;
  return out;
}
std::vector<double> reference_ras_energy_hhe_restricted() {
  std::vector<double> out(2);
  out[0] = -3.26784974;
  out[1] = -2.74836078;
  return out;
}


BOOST_AUTO_TEST_SUITE(TEST_RAS)

BOOST_AUTO_TEST_CASE(FCI_EQUIVALENT) {
    BOOST_CHECK(compare(ras_energy("h2o_sto3g_ras_full"), reference_ras_energy_h2o_full()));
    BOOST_CHECK(compare(ras_energy("hhe_svp_ras_full"), reference_ras_energy_hhe_full()));
}

BOOST_AUTO_TEST_CASE(RESTRICTED) {
    BOOST_CHECK(compare(ras_energy("h2o_sto3g_ras_restricted"), reference_ras_energy_h2o_restricted()));
    BOOST_CHECK(compare(ras_energy("hhe_svp_ras_restricted"), reference_ras_energy_hhe_restricted()));
}

BOOST_AUTO_TEST_SUITE_END()

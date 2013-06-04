//
// BAGEL - Parallel electron correlation program.
// Filename: test_fci.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/fci/harrison.h>
#include <src/fci/knowles.h>

std::vector<double> fci_energy(std::string inp) {

  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = "../../test/" + inp + ".in";
  boost::property_tree::ptree idata;
  boost::property_tree::json_parser::read_json(filename, idata);
  auto keys = idata.get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<Reference> ref;

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    std::string method = iter->second.get<std::string>("title", "");
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(iter->second);

    } else if (method == "hf") {
      auto scf = std::make_shared<SCF>(iter->second, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "rohf"){
      auto scf = std::make_shared<ROHF>(iter->second, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "fci") {
      std::shared_ptr<FCI> fci;
      std::string algorithm = iter->second.get<std::string>("algorithm", "");
      if (algorithm == "harrison") fci = std::make_shared<HarrisonZarrabian>(iter->second, ref);
      else if (algorithm == "knowles") fci = std::make_shared<KnowlesHandy>(iter->second, ref);
      else assert(false);
     
      fci->compute();
      std::cout.rdbuf(backup_stream);
      return fci->energy();
    }
  }
  assert(false);
  return std::vector<double>();
}

std::vector<double> reference_fci_energy() {
  std::vector<double> out(2);
  out[0] = -98.56280393;
  out[1] = -98.36567235;
  return out;
}
std::vector<double> reference_fci_energy2() {
  std::vector<double> out(2);
  out[0] = -3.30315793;
  out[1] = -2.78407879;
  return out;
}


BOOST_AUTO_TEST_SUITE(TEST_NEWFCI)

BOOST_AUTO_TEST_CASE(NEWFCI_2STATE) {
    BOOST_CHECK(compare(fci_energy("hf_sto3g_fci_hz"), reference_fci_energy()));
    BOOST_CHECK(compare(fci_energy("hf_sto3g_fci_kh"), reference_fci_energy()));
    BOOST_CHECK(compare(fci_energy("hhe_svp_fci_hz_trip"), reference_fci_energy2()));
    BOOST_CHECK(compare(fci_energy("hhe_svp_fci_kh_trip"), reference_fci_energy2()));
}

BOOST_AUTO_TEST_SUITE_END()

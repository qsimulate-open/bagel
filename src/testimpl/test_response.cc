//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_response.cc
// Copyright (C) 2017 Toru Shiozaki
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

#include <src/response/cis.h>
#include <src/scf/hf/rhf.h>

using namespace bagel;

std::vector<double> cis_energy(std::string filename, std::string extension = ".json") {
  auto ofs = std::make_shared<std::ofstream>(filename + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << location__ << filename << extension;
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<const Reference> ref;

  std::vector<double> cis_energy;
  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = geom ? std::make_shared<Geometry>(*geom, itree) : std::make_shared<Geometry>(itree);
      if (ref) ref = ref->project_coeff(geom);
    } else if (method == "hf") {
      auto scf = std::make_shared<RHF>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "cis") {
      auto cis = std::make_shared<CIS>(itree, geom, ref);
      cis->compute();
      cis_energy = cis->excitation_energy();
    }
  }
  std::cout.rdbuf(backup_stream);
  return cis_energy;
}

static std::vector<double> hf_svp_cis_ref() {
  return std::vector<double>{0.26095379, 0.26095379, 0.44027041};
}


BOOST_AUTO_TEST_SUITE(TEST_RESPONSE)

BOOST_AUTO_TEST_CASE(CIS) {
    BOOST_CHECK(compare<std::vector<double>>(cis_energy("hf_svp_cis"),  hf_svp_cis_ref(), 1.0e-6));
}

BOOST_AUTO_TEST_SUITE_END()

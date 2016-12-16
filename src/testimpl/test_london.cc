//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_london.cc
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

#include <sstream>
#include <src/scf/giaohf/rhf_london.h>
#include <src/scf/dhf/dirac.h>
#include <src/wfn/relreference.h>

using namespace bagel;

double london_energy(std::string filename) {
  auto ofs = std::make_shared<std::ofstream>(filename + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << location__ << filename << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  std::shared_ptr<const Reference> ref;
  double energy = 0.0;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = geom ? std::make_shared<Geometry>(*geom, itree) : std::make_shared<Geometry>(itree);
      if (ref) ref = ref->project_coeff(geom);
    } else if (method == "hf") {
      auto scf = std::make_shared<RHF_London>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();
      energy = ref->energy(0);
    } else if (method == "dhf") {
      auto rel = std::make_shared<Dirac>(itree, geom, ref);
      rel->compute();
      ref = rel->conv_to_ref();
      energy = ref->energy(0);
    }
  }
  std::cout.rdbuf(backup_stream);
  return energy;
}

BOOST_AUTO_TEST_SUITE(TEST_LONDON)

BOOST_AUTO_TEST_CASE(LONDON) {
  BOOST_CHECK(compare(london_energy("hf_svp_london_hf"),         -99.70397733));
  BOOST_CHECK(compare(london_energy("hf_svp_london_dfhf"),       -99.70391005));
  BOOST_CHECK(compare(london_energy("hf_svp_common_dfhf"),       -98.55057515));
  BOOST_CHECK(compare(london_energy("hf_svp_london_coulomb"),    -99.82461004));
  BOOST_CHECK(compare(london_energy("hf_svp_london_gaunt"),      -99.81352977));
  BOOST_CHECK(compare(london_energy("hcl_svp_london_coulomb"),  -459.16442992));
  BOOST_CHECK(compare(london_energy("hcl_svp_common_coulomb"),  -457.89462676));
  BOOST_CHECK(compare(london_energy("hcl_svp_common_gaunt"),    -454.55922514));
}

BOOST_AUTO_TEST_SUITE_END()

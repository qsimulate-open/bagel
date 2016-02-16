//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_scf.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <sstream>
#include <src/wfn/reference.h>
#include <src/wfn/localization.h>

using namespace bagel;

double localization(std::string filename) {
  auto ofs = std::make_shared<std::ofstream>(filename + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << location__ << filename << ".json";
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

    } else if (method == "localize") {
      if (ref == nullptr) throw std::runtime_error("Localize needs a reference");

      std::string localizemethod = itree->get<std::string>("algorithm", "pm");
      std::shared_ptr<OrbitalLocalization> localization;
      if (localizemethod == "region") {
        localization = std::make_shared<RegionLocalization>(itree, ref);
      }
      else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
        localization = std::make_shared<PMLocalization>(itree, ref);
      else throw std::runtime_error("Unrecognized orbital localization method");

      localization->localize();

      std::cout.rdbuf(backup_stream);
      return localization->metric();
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_LOCALIZE)

BOOST_AUTO_TEST_CASE(PML) {
    BOOST_CHECK(compare(localization("benzene_sto3g_pml"),0.7951349703, 0.000001));
    BOOST_CHECK(compare(localization("watertrimer_sto3g_pml_region"),0.9999109690, 0.000001));
}

BOOST_AUTO_TEST_CASE(REGION) {
    BOOST_CHECK(compare(localization("watertrimer_sto3g_rl"), 0.0));
}

BOOST_AUTO_TEST_SUITE_END()

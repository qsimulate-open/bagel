//
// BAGEL - Parallel electron correlation program.
// Filename: test_scf.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <sstream>
#include <src/wfn/reference.h>
#include <src/molecule/localization.h>

using namespace bagel;

double localization(std::string filename) {
  auto ofs = std::make_shared<std::ofstream>(filename + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
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
    BOOST_CHECK(compare(localization("benzene_sto3g_pml"),13.2770320419, 0.01));
}

BOOST_AUTO_TEST_CASE(REGION) {
    BOOST_CHECK(compare(localization("watertrimer_sto3g_rl"), 0.0));
}

BOOST_AUTO_TEST_SUITE_END()

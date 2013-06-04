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

#include <sstream>
#include <src/wfn/reference.h>
#include <src/util/localization.h>
#include <src/util/lexical_cast.h>

using namespace bagel;

double pm_localization(std::string filename) {
  auto ofs = std::make_shared<std::ofstream>(filename + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".in";
  boost::property_tree::ptree idata;
  boost::property_tree::json_parser::read_json(ss.str(), idata);
  std::shared_ptr<const Reference> ref;
  auto keys = idata.get_child("bagel");

  std::shared_ptr<Geometry> geom;

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    std::string method = iter->second.get<std::string>("title", "");
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(iter->second);

    } else if (method == "hf") {
      auto scf = std::make_shared<SCF>(iter->second, geom);
      scf->compute();
      ref = scf->conv_to_ref();

    } else if (method == "localize") {
      if (ref == nullptr) throw std::runtime_error("Localize needs a reference");

      std::string localizemethod = iter->second.get<std::string>("algorithm", "pm");
      std::shared_ptr<OrbitalLocalization> localization;
      if (localizemethod == "region") {
#if 0
        std::vector<int> sizes;
        auto bound = iter->second.equal_range("region");
        for (auto isizes = bound.first; isizes != bound.second; ++isizes) sizes.push_back(lexical_cast<int>(isizes->second));

        localization = std::make_shared<RegionLocalization>(ref, sizes);
#else
throw std::logic_error("region localization test broken");
#endif
      }   
      else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
        localization = std::make_shared<PMLocalization>(ref);
      else throw std::runtime_error("Unrecognized orbital localization method");

      const int max_iter = iter->second.get<int>("max_iter", 50);
      const double thresh = iter->second.get<double>("thresh", 1.0e-6);

      localization->localize(max_iter, thresh);

      std::cout.rdbuf(backup_stream);
      return localization->metric();
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_LOCALIZE)

BOOST_AUTO_TEST_CASE(PML) {
    BOOST_CHECK(compare(pm_localization("benzene_sto3g_pml"),13.2769000957, 1.0e-5));
}

BOOST_AUTO_TEST_SUITE_END()

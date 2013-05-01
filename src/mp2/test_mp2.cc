//
// BAGEL - Parallel electron correlation program.
// Filename: test_mp2.cc
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


#include <memory>
#include <src/mp2/mp2.h>

double mp2_energy() {

  std::shared_ptr<std::ofstream> ofs(new std::ofstream("benzene_svp_mp2.testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = "../../test/benzene_svp_mp2.in";
  boost::property_tree::ptree idata;
  boost::property_tree::json_parser::read_json(filename, idata);
  auto keys = idata.get_child("bagel");
  std::shared_ptr<Geometry> geom;

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    std::string method = iter->second.get<std::string>("title", "");
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);

    if (method == "molecule") {
      geom = std::shared_ptr<Geometry>(new Geometry(iter->second));

    } else if (method == "mp2") {
      std::shared_ptr<MP2> mp2(new MP2(iter->second, geom));
      mp2->compute();

      std::cout.rdbuf(backup_stream);
      return mp2->energy();
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_MP2)

BOOST_AUTO_TEST_CASE(MP2) {
    BOOST_CHECK(compare(mp2_energy(), -231.31440958));
}

BOOST_AUTO_TEST_SUITE_END()

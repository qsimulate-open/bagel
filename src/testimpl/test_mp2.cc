//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_mp2.cc
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


#include <memory>
#include <src/pt2/mp2/mp2.h>

double mp2_energy(const std::string& job) {

  auto ofs = std::make_shared<std::ofstream>(job + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = location__ + job + ".json";
  auto idata = std::make_shared<const PTree>(filename);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "mp2") {
      auto mp2 = std::make_shared<MP2>(itree, geom);
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
    BOOST_CHECK(compare(mp2_energy("benzene_svp_mp2"),      -231.31440958));
    BOOST_CHECK(compare(mp2_energy("benzene_svp_mp2_aux"),  -231.31450878));
}

BOOST_AUTO_TEST_SUITE_END()

//
// BAGEL - Parallel electron correlation program.
// Filename: test_nevpt2.cc
// Copyright (C) 2014 Toru Shiozaki
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


#include <memory>
#include <src/nevpt2/nevpt2.h>

double nevpt2_energy(const std::string job) {

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

    } else if (method == "nevpt2") {
      auto nevpt2 = std::make_shared<NEVPT2>(itree, geom);
      nevpt2->compute();

      std::cout.rdbuf(backup_stream);
      return nevpt2->energy();
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_NEVPT2)

BOOST_AUTO_TEST_CASE(NEVPT2) {
    BOOST_CHECK(compare(nevpt2_energy("h2o_svp_nevpt2"),      -76.0205249495));
}

BOOST_AUTO_TEST_SUITE_END()

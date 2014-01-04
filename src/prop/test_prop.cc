//
// BAGEL - Parallel electron correlation program.
// Filename: test_prop.cc
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
#include <src/prop/multipole.h>
#include <src/scf/scf.h>
#include <src/scf/rohf.h>
#include <src/scf/uhf.h>
#include <src/wfn/reference.h>

std::vector<double> multipole(std::string filename) {
  auto ofs = std::make_shared<std::ofstream>(filename + "_multipole.testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<SCF>(itree, geom);
      scf->compute();
      std::shared_ptr<const Matrix> dtot = scf->coeff()->form_density_rhf(scf->nocc());

      Multipole multipole(geom, dtot, 2);
      std::vector<double> d = multipole.compute();
      std::cout.rdbuf(backup_stream);
      return d;
    }
  }
  assert(false);
  return std::vector<double>();
}

static std::vector<double> hf_svp_dfhf_multipole_ref() {
  return std::vector<double>{0.0, 0.0, 1.055510, -4.236243, 0.000000, -4.236243, -0.000000, -0.000000, -1.532119};
}

BOOST_AUTO_TEST_SUITE(TEST_PROP)

BOOST_AUTO_TEST_CASE(MULTIPOLE) {
    BOOST_CHECK(compare<std::vector<double>>(multipole("hf_svp_dfhf"),        hf_svp_dfhf_multipole_ref(), 1.0e-6));
}

BOOST_AUTO_TEST_SUITE_END()

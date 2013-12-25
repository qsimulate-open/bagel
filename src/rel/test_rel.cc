//
// BAGEL - Parallel electron correlation program.
// Filename: test_rel.cc
// Copyright (C) 2013 Toru Shiozaki
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
#include <src/rel/dirac.h>
#include <src/wfn/reference.h>

using namespace bagel;

double rel_energy(std::string filename) {
  auto ofs = std::make_shared<std::ofstream>(filename + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  std::shared_ptr<const Reference> ref_;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<SCF>(itree, geom);
      scf->compute();
      ref_ = scf->conv_to_ref();
    } else if (method == "dhf") {
      auto rel = std::make_shared<Dirac>(itree, geom, ref_);
      rel->compute();
      ref_ = rel->conv_to_ref();
      std::cout.rdbuf(backup_stream);
      return ref_->energy();
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_REL)

BOOST_AUTO_TEST_CASE(DIRAC_FOCK) {
    BOOST_CHECK(compare(rel_energy("hf_svp_coulomb"),        -99.93791152));
    BOOST_CHECK(compare(rel_energy("hf_svp_gaunt"),          -99.92699858));
    BOOST_CHECK(compare(rel_energy("hf_svp_breit"),          -99.92755305));
}

BOOST_AUTO_TEST_SUITE_END()

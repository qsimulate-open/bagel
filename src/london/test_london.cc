//
// BAGEL - Parallel electron correlation program.
// Filename: test_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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
#include <src/london/scf_london.h>
#include <src/london/dirac_london.h>
#include <src/london/reference_london.h>

using namespace bagel;

double london_energy(std::string filename) {
  auto ofs = std::make_shared<std::ofstream>(filename + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry_London> cgeom;

  std::shared_ptr<const Reference> ref_;
  double energy = 0.0;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      cgeom = std::make_shared<Geometry_London>(itree);
    } else if (method == "hf") {
      auto scf = std::make_shared<SCF_London>(itree, cgeom, ref_);
      scf->compute();
      ref_ = scf->conv_to_ref();
      energy = ref_->energy();
    } else if (method == "dhf") {
      auto rel = std::make_shared<Dirac_London>(itree, cgeom, ref_);
      rel->compute();
      ref_ = rel->conv_to_ref();
      //return ref_->energy(); //TODO Set up reference
      energy =  rel->energy();
    }
  }
  std::cout.rdbuf(backup_stream);
  return energy;
}

BOOST_AUTO_TEST_SUITE(TEST_LONDON)

BOOST_AUTO_TEST_CASE(LONDON) {
  BOOST_CHECK(compare(london_energy("hf_svp_london_hf"),      -99.70397733));
  BOOST_CHECK(compare(london_energy("hf_svp_london_dfhf"),    -99.70391005));
  BOOST_CHECK(compare(london_energy("hf_svp_london_coulomb"), -99.82459461));
}

BOOST_AUTO_TEST_SUITE_END()

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
#include <src/rel/dirac.h>
#include <src/wfn/reference.h>

using namespace bagel;

double rel_energy(std::string filename) {
  std::shared_ptr<std::ofstream> ofs(new std::ofstream(filename + ".testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".in";
  std::shared_ptr<InputData> idata(new InputData(ss.str()));
  std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
  std::list<std::pair<std::string, std::multimap<std::string, std::string>>> keys = idata->data();

  std::shared_ptr<Reference> ref_;

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "df-hf") {
      std::shared_ptr<SCF<1>> scf(new SCF<1>(iter->second, geom));
      scf->compute();
      ref_ = scf->conv_to_ref();
    } else if (iter->first == "dhf") {
      std::shared_ptr<Dirac> rel(new Dirac(iter->second, geom, ref_));
      rel->compute();
      std::shared_ptr<RelReference> ref = rel->conv_to_ref();
      std::cout.rdbuf(backup_stream);
      return ref->energy();
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_REL)

BOOST_AUTO_TEST_CASE(BREIT) {
    BOOST_CHECK(compare(rel_energy("hf_svp_breit"),          -99.92755305));
}

BOOST_AUTO_TEST_SUITE_END()

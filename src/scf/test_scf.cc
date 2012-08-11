//
// Newint - Parallel electron correlation program.
// Filename: test_scf.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <sstream>
#include <src/scf/scf.h>
#include <src/scf/rohf.h>
#include <src/scf/uhf.h>
#include <src/wfn/reference.h>

double scf_energy(std::string filename) {
  std::shared_ptr<std::ofstream> ofs(new std::ofstream(filename + ".testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".in";
  std::shared_ptr<InputData> idata(new InputData(ss.str()));
  stack = new StackMem(static_cast<size_t>(1000000LU));
  std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
  std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "df-hf") {
      std::shared_ptr<SCF<1> > scf(new SCF<1>(iter->second, geom));
      scf->compute();
      std::shared_ptr<Reference> ref = scf->conv_to_ref();

      delete stack;
      std::cout.rdbuf(backup_stream);
      return ref->energy();
    } else if (iter->first == "df-uhf") {
      std::shared_ptr<UHF> scf(new UHF(iter->second, geom));
      scf->compute();
      std::shared_ptr<Reference> ref = scf->conv_to_ref();

      delete stack;
      std::cout.rdbuf(backup_stream);
      return ref->energy();
    } else if (iter->first == "df-rohf") {
      std::shared_ptr<ROHF> scf(new ROHF(iter->second, geom));
      scf->compute();
      std::shared_ptr<Reference> ref = scf->conv_to_ref();

      delete stack;
      std::cout.rdbuf(backup_stream);
      return ref->energy();
    } else if (iter->first == "hf") {
      std::shared_ptr<SCF<0> > scf(new SCF<0>(iter->second, geom));
      scf->compute();
      std::shared_ptr<Reference> ref = scf->conv_to_ref();

      delete stack;
      std::cout.rdbuf(backup_stream);
      return ref->energy();
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_SCF)
 
BOOST_AUTO_TEST_CASE(DF_HF) {
    BOOST_CHECK(compare(scf_energy("hf_svp_hf"),          -99.84779026));
    BOOST_CHECK(compare(scf_energy("hf_svp_dfhf"),        -99.84772354));
    BOOST_CHECK(compare(scf_energy("hf_svp_dfhf_ext"),    -99.83765614));
    BOOST_CHECK(compare(scf_energy("hf_svp_dfhf_cart"),   -99.84911270));
    BOOST_CHECK(compare(scf_energy("hf_svp_dfhf_charge"), -99.77381455));
    BOOST_CHECK(compare(scf_energy("oh_svp_uhf"),         -75.28410147));
    BOOST_CHECK(compare(scf_energy("hc_svp_rohf"),        -38.16810629));
}
 
BOOST_AUTO_TEST_SUITE_END()

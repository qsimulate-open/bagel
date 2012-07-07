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
#include <src/wfn/reference.h>

double scf_energy(std::string filename) {
  std::shared_ptr<std::ofstream> ofs(new std::ofstream(filename + ".testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".in";
  std::shared_ptr<InputData> idata(new InputData(ss.str()));
  stack = new StackMem(static_cast<size_t>(1000000LU));
  std::shared_ptr<Geometry> geom(new Geometry(idata));
  std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "df-hf") {
      std::shared_ptr<SCF<1> > scf(new SCF<1>(iter->second, geom));
      scf->compute();
      std::shared_ptr<Reference> ref = scf->conv_to_ref();

      delete stack;
      std::cout.rdbuf(backup_stream);
      return ref->energy();
    }
  }
  assert(false);
}

BOOST_AUTO_TEST_SUITE(TEST_SCF)
 
BOOST_AUTO_TEST_CASE(DF_HF) {
    BOOST_CHECK(compare(scf_energy("hf_svp_dfhf"),     -98.49666065));
    BOOST_CHECK(compare(scf_energy("hf_svp_dfhf_ext"), -98.49306775));
}
 
BOOST_AUTO_TEST_SUITE_END()

//
// Newint - Parallel electron correlation program.
// Filename: test_casscf.cc
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
#include <src/casscf/superci.h>
#include <src/wfn/reference.h>

double cas_energy(std::string filename) {
  std::shared_ptr<std::ofstream> ofs(new std::ofstream(filename + ".testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".in";
  std::shared_ptr<InputData> idata(new InputData(ss.str()));
  stack = new StackMem(static_cast<size_t>(1000000LU));
  std::shared_ptr<Geometry> geom(new Geometry(idata));
  std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "casscf") {
      std::string algorithm = read_input<std::string>(iter->second, "algorithm", ""); 
      if (algorithm == "superci" || algorithm == "") {
        std::shared_ptr<CASSCF> cas(new SuperCI(iter->second, geom));
        cas->compute();
        std::shared_ptr<const Reference> ref = cas->conv_to_ref();

        delete stack;
        std::cout.rdbuf(backup_stream);
        return ref->energy();
      }
    }
  }
  assert(false);
}

BOOST_AUTO_TEST_SUITE(TEST_CASSCF)
 
BOOST_AUTO_TEST_CASE(DF_CASSCF) {
    BOOST_CHECK(compare(cas_energy("lif_svp_cas22"),      -106.70563743));
}
 
BOOST_AUTO_TEST_SUITE_END()

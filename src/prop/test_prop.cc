//
// Newint - Parallel electron correlation program.
// Filename: test_prop.cc
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

std::array<double,3> dipole(std::string filename) {
  std::shared_ptr<std::ofstream> ofs(new std::ofstream(filename + "_dipole.testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << filename << ".in";
  std::shared_ptr<InputData> idata(new InputData(ss.str()));
  stack = new StackMem();
  std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
  std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "df-hf") {
      std::shared_ptr<SCF<1> > scf(new SCF<1>(iter->second, geom));
      scf->compute();
      std::shared_ptr<const Matrix1e> dtot = scf->coeff()->form_density_rhf(scf->nocc());

      Dipole dipole(geom, dtot);
      std::array<double,3> d = dipole.compute();
      delete stack;
      std::cout.rdbuf(backup_stream);
      return d;
    }
  }
  assert(false);
  return std::array<double,3>();
}

static std::array<double,3> hf_svp_dfhf_dipole_ref() {
  return std::array<double,3>{{0.0, 0.0, 1.055510}};
}

typedef std::array<double,3> ARRAY;

BOOST_AUTO_TEST_SUITE(TEST_PROP)
 
BOOST_AUTO_TEST_CASE(DIPOLE) {
    BOOST_CHECK(compare<ARRAY>(dipole("hf_svp_dfhf"),        hf_svp_dfhf_dipole_ref(), 1.0e-6));
}
 
BOOST_AUTO_TEST_SUITE_END()

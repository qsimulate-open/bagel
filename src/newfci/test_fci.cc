//
// Newint - Parallel electron correlation program.
// Filename: test_fci.cc
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


#include <src/newfci/fci.h>

std::vector<double> fci_energy() {

  std::shared_ptr<std::ofstream> ofs(new std::ofstream("hf_sto3g_fci2.testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::shared_ptr<InputData> idata(new InputData("../../test/hf_sto3g_fci2.in"));
  stack = new StackMem(static_cast<size_t>(1000000LU));
  std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
  std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

  std::shared_ptr<Reference> ref;

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "df-hf") {
      std::shared_ptr<SCF<1> > scf(new SCF<1>(iter->second, geom));
      scf->compute();
      ref = scf->conv_to_ref();

    }
    if (iter->first == "fci") {
      std::shared_ptr<NewFCI> fci(new NewFCI(iter->second, ref));
      fci->compute();

      delete stack;
      std::cout.rdbuf(backup_stream);
      return fci->energy();
    }
  }
  assert(false);
  return std::vector<double>();
}

std::vector<double> reference_fci_energy() {
  std::vector<double> out(2);
  out[0] = -98.56280393;
  out[1] = -98.36567235;
  return out;
}
 
BOOST_AUTO_TEST_SUITE(TEST_NewFCI)
 
BOOST_AUTO_TEST_CASE(NewFCI_2STATE) {
    BOOST_CHECK(compare(fci_energy(), reference_fci_energy()));
}
 
BOOST_AUTO_TEST_SUITE_END()

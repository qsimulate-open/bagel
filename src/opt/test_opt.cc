//
// Newint - Parallel electron correlation program.
// Filename: test_opt.cc
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

#include <src/opt/opt.h>
#include <src/scf/scf.h>
#include <src/wfn/reference.h>

std::vector<double> scf_opt() {

  std::shared_ptr<std::ofstream> ofs(new std::ofstream("hf_svp_dfhf_opt.testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::shared_ptr<InputData> idata(new InputData("../../test/hf_svp_dfhf_opt.in"));
  stack = new StackMem(static_cast<size_t>(1000000LU));
  std::shared_ptr<Geometry> geom(new Geometry(idata));
  std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "df-hf-opt") {
      std::shared_ptr<Opt<SCF<1> > > opt(new Opt<SCF<1> >(idata, iter->second, geom));
      for (int i = 0; i != 20; ++i)
        if (opt->next()) break;

      delete stack;
      std::cout.rdbuf(backup_stream);
      return opt->geometry()->xyz();
    }
  }
  assert(false);
}
std::vector<double> reference_scf_opt() {
  std::vector<double> out(6);
  out[2] = 1.864207;
  out[5] = 0.162365;
  return out;
}

std::vector<double> mp2_opt() {

  std::shared_ptr<std::ofstream> ofs(new std::ofstream("hf_svp_mp2_opt.testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::shared_ptr<InputData> idata(new InputData("../../test/hf_svp_mp2_opt.in"));
  stack = new StackMem(static_cast<size_t>(1000000LU));
  std::shared_ptr<Geometry> geom(new Geometry(idata));
  std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "mp2-opt") {
      std::shared_ptr<Opt<MP2Grad> > opt(new Opt<MP2Grad>(idata, iter->second, geom));
      for (int i = 0; i != 20; ++i)
        if (opt->next()) break;

      delete stack;
      std::cout.rdbuf(backup_stream);
      return opt->geometry()->xyz();
    }
  }
  assert(false);
}
std::vector<double> reference_mp2_opt() {
  std::vector<double> out(6);
  out[2] = 1.981406;
  out[5] = 0.245166;
  return out;
}
 
BOOST_AUTO_TEST_SUITE(TEST_OPT)
 
BOOST_AUTO_TEST_CASE(DF_HF_Opt) {
    BOOST_CHECK(compare(scf_opt(), reference_scf_opt(), 1.0e-6));
}
BOOST_AUTO_TEST_CASE(MP2_Opt) {
    BOOST_CHECK(compare(mp2_opt(), reference_mp2_opt(), 1.0e-6));
}
 
BOOST_AUTO_TEST_SUITE_END()

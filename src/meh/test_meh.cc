//
// BAGEL - Parallel electron correlation program.
// Filename: test_meh.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <src/meh/meh.h>
#include <src/dimer/dimer.h>

double meh_energy(std::string inp) {

  std::shared_ptr<std::ofstream> ofs(new std::ofstream(inp + ".testout", std::ios::trunc));
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::shared_ptr<InputData> idata(new InputData("../../test/" + inp + ".in"));
  std::shared_ptr<Geometry> geom(new Geometry(idata->get_input("molecule")));
  std::list<std::pair<std::string, std::multimap<std::string, std::string>>> keys = idata->data();

  std::shared_ptr<Reference> ref;
  std::shared_ptr<Dimer> dimer;

  for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    if (iter->first == "df-hf") {
      std::shared_ptr<SCF<1>> scf(new SCF<1>(iter->second, geom));
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (iter->first == "dimerize") { // dimerize forms the dimer object, does a scf calculation, and then localizes
      std::multimap<std::string,std::string> dimdata = iter->second;

      std::string form = read_input<std::string>(dimdata, "form", "displace");
      if (form == "d" || form == "disp" || form == "displace") {
        double scale = (read_input<bool>(dimdata,"angstrom",false) ? ang2bohr__ : 1.0 ) ; 

        double dx = read_input<double>(dimdata,"dx",0.0) * scale;
        double dy = read_input<double>(dimdata,"dy",0.0) * scale;
        double dz = read_input<double>(dimdata,"dz",0.0) * scale;
        std::array<double,3> disp = {{dx,dy,dz}};

        if (static_cast<bool>(ref)) {
          dimer = std::make_shared<Dimer>(ref,disp);
        }   
        else {
          throw std::runtime_error("dimerize needs a reference calculation (for now)");
        }   
      }   
      dimer->scf(iter->second);

      *geom = *dimer->sgeom();
      ref = dimer->sref();
    } else if (iter->first == "meh") {
      std::shared_ptr<DimerCISpace> cispace = dimer->compute_cispace(iter->second);

      auto meh = std::make_shared<MultiExcitonHamiltonian>(dimer, cispace);
      meh->compute();

      std::cout.rdbuf(backup_stream);
      return meh->energy(0);
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_MEH)

BOOST_AUTO_TEST_CASE(MEH_GROUND_STATE) {
    BOOST_CHECK(compare(meh_energy("benzene_sto3g_meh"), -459.319548355058));
}

BOOST_AUTO_TEST_SUITE_END()

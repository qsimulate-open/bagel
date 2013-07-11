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

  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << "../../test/" << inp << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  std::shared_ptr<Reference> ref;
  std::shared_ptr<Dimer> dimer;

  std::map<std::string, std::shared_ptr<const void>> saved;

  for (auto& itree : *keys) {
    std::string method = itree->get<std::string>("title", "");
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<SCF>(itree, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "dimerize") { // dimerize forms the dimer object, does a scf calculation, and then localizes

      std::string form = itree->get<std::string>("form", "displace");
      if (form == "d" || form == "disp" || form == "displace") {
        double scale = (itree->get<bool>("angstrom",false) ? ang2bohr__ : 1.0 ) ;

        double dx = itree->get<double>("dx",0.0) * scale;
        double dy = itree->get<double>("dy",0.0) * scale;
        double dz = itree->get<double>("dz",0.0) * scale;
        std::array<double,3> disp = {{dx,dy,dz}};

        if (static_cast<bool>(ref)) {
          dimer = std::make_shared<Dimer>(ref,disp);
        } else {
          throw std::runtime_error("dimerize needs a reference calculation (for now)");
        }
      }
      else if (form == "r" || form == "refs") {
        std::vector<std::shared_ptr<const Reference>> dimer_refs;
        auto units = itree->get_vector<std::string>("refs", 2);
        for (auto& ikey : units) {
          auto tmp = saved.find(ikey);
          if (tmp == saved.end()) throw std::runtime_error(std::string("No reference found with name: ") + ikey);
          else dimer_refs.push_back(std::static_pointer_cast<const Reference>(tmp->second));
        }

        dimer = std::make_shared<Dimer>(dimer_refs.at(0), dimer_refs.at(1));
      }
      dimer->scf(itree);

      *geom = *dimer->sgeom();
      ref = dimer->sref();
    } else if (method == "meh") {
      std::shared_ptr<DimerCISpace> cispace = dimer->compute_cispace(itree);

      auto meh = std::make_shared<MultiExcitonHamiltonian>(itree, dimer, cispace);
      meh->compute();

      std::cout.rdbuf(backup_stream);
      return meh->energy(0);
    }

    std::string saveref = itree->get<std::string>("saveref", "");
    if (saveref != "") { saved.insert(std::make_pair(saveref, std::static_pointer_cast<const void>(ref))); }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_MEH)

BOOST_AUTO_TEST_CASE(MEH_GROUND_STATE) {
    BOOST_CHECK(compare(meh_energy("benzene_sto3g_meh_stack"), -459.319548355058));
    BOOST_CHECK(compare(meh_energy("benzene_sto3g_meh_T"), -459.28040890));
}

BOOST_AUTO_TEST_SUITE_END()

//
// BAGEL - Parallel electron correlation program.
// Filename: test_asd_dmrg.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/multisite/multisite.h>
#include <src/asd/dmrg/rasd.h>


double asd_dmrg_energy(std::string inp) {
  Muffle hide_cout(inp + ".testout");

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << location__ << inp << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  std::shared_ptr<const Reference> ref;
  std::shared_ptr<MultiSite> multisite;

  std::map<std::string, std::shared_ptr<const Reference>> saved;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<RHF>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "multisite") {
      std::vector<std::shared_ptr<const Reference>> site_refs;
      auto sitenames = itree->get_vector<std::string>("refs");
      for (auto& s : sitenames)
        site_refs.push_back(saved.at(s));
      auto ms = std::make_shared<MultiSite>(itree, site_refs);
      ms->scf(itree);
      multisite = ms;
      ref = ms->conv_to_ref();
      *geom = *ref->geom();
    } else if (method == "asd_dmrg") {
        if (!multisite)
          throw std::runtime_error("multisite must be called before asd_dmrg");
        auto asd = std::make_shared<RASD>(itree, multisite);
        asd->compute();
        return asd->energies(0);
    }

    std::string saveref = itree->get<std::string>("saveref", "");
    if (saveref != "") { saved.emplace(saveref, ref); }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_ASD_DMRG)

BOOST_AUTO_TEST_CASE(RAS) {
    BOOST_CHECK(compare(asd_dmrg_energy("he3_svp_rasd-dmrg"), -8.53974980, 1.0e-8));
}

BOOST_AUTO_TEST_SUITE_END()

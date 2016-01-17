//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_asd.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <src/asd/construct_asd.h>

double asd_energy(std::string inp) {
  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << location__ << inp << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  std::shared_ptr<const Reference> ref;
  std::shared_ptr<Dimer> dimer;

  std::map<std::string, std::shared_ptr<const void>> saved;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<RHF>(itree, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "dimerize") { // dimerize forms the dimer object, does a scf calculation, and then localizes
      std::string form = itree->get<std::string>("form", "displace");
      if (form == "d" || form == "disp" || form == "displace") {
        if (static_cast<bool>(ref))
          dimer = std::make_shared<Dimer>(itree, ref);
        else
          throw std::runtime_error("dimerize needs a reference calculation (for now)");
      }
      else if (form == "r" || form == "refs") {
        std::vector<std::shared_ptr<const Reference>> dimer_refs;
        auto units = itree->get_vector<std::string>("refs", 2);
        for (auto& ikey : units) {
          auto tmp = saved.find(ikey);
          if (tmp == saved.end()) throw std::runtime_error(std::string("No reference found with name: ") + ikey);
          else dimer_refs.push_back(std::static_pointer_cast<const Reference>(tmp->second));
        }

        dimer = std::make_shared<Dimer>(itree, dimer_refs.at(0), dimer_refs.at(1));
      }
      dimer->scf(itree);

      *geom = *dimer->sgeom();
      ref = dimer->sref();
    } else if (method == "asd") {
      auto asd = construct_ASD(itree, dimer);
      asd->compute();

      std::cout.rdbuf(backup_stream);
      return asd->energy(0);
    }

    std::string saveref = itree->get<std::string>("saveref", "");
    if (saveref != "") { saved.emplace(saveref, std::static_pointer_cast<const void>(ref)); }
  }
  assert(false);
  return 0.0;
}

std::vector<double> asd_models(std::string inp) {
  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << location__ << inp << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  std::shared_ptr<const Reference> ref;
  std::shared_ptr<Dimer> dimer;

  std::map<std::string, std::shared_ptr<const void>> saved;

  std::vector<double> models;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<RHF>(itree, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "dimerize") { // dimerize forms the dimer object, does a scf calculation, and then localizes
      std::string form = itree->get<std::string>("form", "displace");
      if (form == "d" || form == "disp" || form == "displace") {
        if (static_cast<bool>(ref))
          dimer = std::make_shared<Dimer>(itree, ref);
        else
          throw std::runtime_error("dimerize needs a reference calculation (for now)");
      }
      else if (form == "r" || form == "refs") {
        std::vector<std::shared_ptr<const Reference>> dimer_refs;
        auto units = itree->get_vector<std::string>("refs", 2);
        for (auto& ikey : units) {
          auto tmp = saved.find(ikey);
          if (tmp == saved.end()) throw std::runtime_error(std::string("No reference found with name: ") + ikey);
          else dimer_refs.push_back(std::static_pointer_cast<const Reference>(tmp->second));
        }

        dimer = std::make_shared<Dimer>(itree, dimer_refs.at(0), dimer_refs.at(1));
      }
      dimer->scf(itree);

      *geom = *dimer->sgeom();
      ref = dimer->sref();
    } else if (method == "asd") {
      auto asd = construct_ASD(itree, dimer);
      asd->compute();

      std::shared_ptr<Matrix> dirmodel = asd->model(0).first;
      models.insert(models.end(), dirmodel->data(), dirmodel->data()+dirmodel->size());
      models[1] = std::abs(models[1]); models[2] = std::abs(models[2]);

      std::shared_ptr<Matrix> ptmodel = asd->model(0).second;
      models.insert(models.end(), ptmodel->data(), ptmodel->data()+ptmodel->size());
      models[5] = std::abs(models[5]); models[6] = std::abs(models[6]);
    }

    std::string saveref = itree->get<std::string>("saveref", "");
    if (saveref != "") { saved.emplace(saveref, std::static_pointer_cast<const void>(ref)); }
  }

  std::cout.rdbuf(backup_stream);
  return models;
}

BOOST_AUTO_TEST_SUITE(TEST_ASD)

BOOST_AUTO_TEST_CASE(CAS) {
    BOOST_CHECK(compare(asd_energy("benzene_sto3g_asd_stack"), -459.40037137, 1.0e-6));
    BOOST_CHECK(compare(asd_energy("benzene_sto3g_asd_T"), -459.36294726, 1.0e-6));
}

BOOST_AUTO_TEST_CASE(RAS) {
    BOOST_CHECK(compare(asd_models("benzene_sto3g_asd-ras_stack"),
      std::vector<double>{{ 0.000000000, 0.000047060, 0.000047060,  0.000000000,
                           -0.001758099, 0.001599564, 0.001599564, -0.001758099 }}, 1.0e-6));
}

BOOST_AUTO_TEST_SUITE_END()

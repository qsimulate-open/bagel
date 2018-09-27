//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_zfs.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/multi/zcasscf/zcassecond.h>
#include <src/prop/pseudospin/pseudospin.h>

// Checking g-values to the 8th decimal place is not ideal
const double gscale = 1.0e-4;

std::vector<double> zfs_param(std::string inp) {

  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());
  std::vector<double> out;

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = location__ + inp + ".json";
  auto idata = std::make_shared<const PTree>(filename);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<const Reference> ref;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<RHF>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();

    } else if (method == "dhf") {
      auto scf = std::make_shared<Dirac>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();

    } else if (method == "zcasscf") {
      std::string algorithm = itree->get<std::string>("algorithm", "");
      if (algorithm == "second" || algorithm == "") {
        auto zcas = std::make_shared<ZCASSecond>(itree, geom, ref);
        zcas->compute();
        ref = zcas->conv_to_ref();
        std::shared_ptr<const PTree> aniso_data = itree->get_child_optional("aniso");
        if (aniso_data) {
          const int nspin = aniso_data->get<int>("nspin", 3);
          auto ps = std::make_shared<Pseudospin>(nspin, geom->relativistic(/*gaunt*/false), zcas->fci()->conv_to_ciwfn(), aniso_data);
          ps->compute(zcas->energy(), zcas->coeff()->active_part());
          out.clear();
          out.push_back(ps->gval(0) * gscale);
          out.push_back(ps->gval(1) * gscale);
          out.push_back(ps->gval(2) * gscale);
          out.push_back(ps->dval());
          out.push_back(ps->eval());
        }
      }
    }
  }
  std::cout.rdbuf(backup_stream);
  return out;
}

std::vector<double> reference_pseudospin_parameters() {
  std::vector<double> out(5);
  out[0] =  2.00430269 * gscale;
  out[1] =  2.00430219 * gscale;
  out[2] =  2.00221459 * gscale;
  out[3] =  0.00000164;
  out[4] =  0.00000000;
  return out;
}

BOOST_AUTO_TEST_SUITE(TEST_PSEUDOSPIN)

BOOST_AUTO_TEST_CASE(ZFSPLITTING) {
  BOOST_CHECK(compare(zfs_param("o2_321g_zcasscf_pseudospin"), reference_pseudospin_parameters()));
}

BOOST_AUTO_TEST_SUITE_END()

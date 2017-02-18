//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_casscf.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/multi/casscf/superci.h>
#include <src/multi/casscf/casbfgs.h>
#include <src/multi/casscf/cashybrid.h>
#include <src/multi/casscf/cassecond.h>

double cas_energy(std::string filename) {
  auto ofs = std::make_shared<std::ofstream>(filename + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << location__ << filename << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<const Reference> ref;

  double energy = 0.0;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);
    } else if (method == "hf") {
      auto hf = std::make_shared<RHF>(itree, geom);
      hf->compute();
      ref = hf->conv_to_ref();

    } else if (method == "casscf") {
      std::string algorithm = itree->get<std::string>("algorithm", "");
      if (algorithm == "superci") {
        auto cas = std::make_shared<SuperCI>(itree, geom, ref);
        cas->compute();
        ref = cas->conv_to_ref();
        energy = ref->energy(0);
      } else if (algorithm == "second") {
        auto cas = std::make_shared<CASSecond>(itree, geom, ref);
        cas->compute();
        ref = cas->conv_to_ref();
        energy = ref->energy(0);
      } else if (algorithm == "hybrid") {
        auto cas = std::make_shared<CASHybrid>(itree, geom, ref);
        cas->compute();
        ref = cas->conv_to_ref();
        energy = ref->energy(0);
      } else if (algorithm == "bfgs") {
        auto cas = std::make_shared<CASBFGS>(itree, geom, ref);
        cas->compute();
        ref = cas->conv_to_ref();
        energy = ref->energy(0);
      }
    }
  }
  assert(energy != 0.0);
  std::cout.rdbuf(backup_stream);
  return energy;
}

BOOST_AUTO_TEST_SUITE(TEST_CASSCF)

BOOST_AUTO_TEST_CASE(DF_CASSCF) {
    BOOST_CHECK(compare(cas_energy("h2o_svp_cassecond"),    -76.00368392));
    BOOST_CHECK(compare(cas_energy("lif_svp_cas22"),        -106.70563743));
    BOOST_CHECK(compare(cas_energy("li2_tzvpp_cas43"),      -14.87300366));
    BOOST_CHECK(compare(cas_energy("lih_tzvpp_cas22"),      -7.98191070));
//    BOOST_CHECK(compare(cas_energy("crco6_sto3g_cas66"),    -1699.66508838));
}

BOOST_AUTO_TEST_SUITE_END()

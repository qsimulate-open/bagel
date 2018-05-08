//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_zcasscf.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson E. Bates  <jefferson.bates@northwestern.edu>
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

double relcas_energy(std::string inp) {

  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = location__ + inp + ".json";
  auto idata = std::make_shared<const PTree>(filename);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<const Reference> ref;
  double energy = 0.0;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = geom ? std::make_shared<Geometry>(*geom, itree) : std::make_shared<Geometry>(itree);
      if (itree->get<bool>("restart", false))
        ref.reset();
      if (ref) ref = ref->project_coeff(geom);

    } else if (method == "hf") {
      auto scf = std::make_shared<RHF>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();

    } else if (method == "dhf") {
      auto scf = std::make_shared<Dirac>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();

    } else if (method == "casscf") {
      auto cas = std::make_shared<CASSecond>(itree, geom, ref);
      cas->compute();
      ref = cas->conv_to_ref();
      energy = ref->energy(0);

    } else if (method == "zcasscf") {
      std::string algorithm = itree->get<std::string>("algorithm", "");
      if (algorithm == "second" || algorithm == "") {
        auto zcas = std::make_shared<ZCASSecond>(itree, geom, ref);
        zcas->compute();
        ref = zcas->conv_to_ref();
        energy = ref->energy(ref->nstate()-1);
      } else {
        throw std::logic_error("unknown algorithm");
      }

#ifndef DISABLE_SERIALIZATION
    } else if (method == "save_ref") {
      const std::string name = itree->get<std::string>("file", "");
      assert (name != "");
      OArchive archive(name);
      archive << ref;

    } else if (method == "load_ref") {
      const std::string name = itree->get<std::string>("file", "");
      assert (name != "");
      IArchive archive(name);
      std::shared_ptr<Reference> ptr;
      archive >> ptr;
      ref = std::shared_ptr<Reference>(ptr);
#endif
    }
  }
  std::cout.rdbuf(backup_stream);
  return energy;
}

BOOST_AUTO_TEST_SUITE(TEST_RELCAS)

BOOST_AUTO_TEST_CASE(ZCASSCF) {
  BOOST_CHECK(compare(relcas_energy("hf_svp_second_coulomb"),   -99.90025083));
#ifndef HAVE_MPI_H
  BOOST_CHECK(compare(relcas_energy("he_tzvpp_second_coulomb"),    -2.875647885));
#endif
  BOOST_CHECK(compare(relcas_energy("nh_svp_triplet_gaunt"),    -54.93378556));
  BOOST_CHECK(compare(relcas_energy("o2_svp_triplet_breit"),   -149.56647946));
#ifndef DISABLE_SERIALIZATION
  BOOST_CHECK(compare(relcas_energy("hf_tzvpp_zcasscf_saveref"), -100.03016820));
#endif
}

BOOST_AUTO_TEST_SUITE_END()

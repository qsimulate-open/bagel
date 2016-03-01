//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_fci.cc
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


#include <src/ci/fci/harrison.h>
#include <src/ci/fci/knowles.h>
#include <src/ci/fci/distfci.h>

std::vector<double> fci_energy(std::string inp) {

  auto ofs = std::make_shared<std::ofstream>(inp + ".testout", std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  // a bit ugly to hardwire an input file, but anyway...
  std::string filename = location__ + inp + ".json";
  auto idata = std::make_shared<const PTree>(filename);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;
  std::shared_ptr<const Reference> ref;

  std::vector<double> result;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<RHF>(itree, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "rohf"){
      auto scf = std::make_shared<ROHF>(itree, geom);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "fci") {
      std::shared_ptr<FCI> fci;
      std::shared_ptr<DistFCI> dfci;
      std::string algorithm = itree->get<std::string>("algorithm", "knowles");
      if (algorithm == "harrison") fci = std::make_shared<HarrisonZarrabian>(itree, geom, ref);
      else if (algorithm == "knowles") fci = std::make_shared<KnowlesHandy>(itree, geom, ref);
      else if (algorithm == "dist" || algorithm == "parallel")
        dfci = std::make_shared<DistFCI>(itree, geom, ref);
      else assert(false);

      if (fci) fci->compute();
      else if (dfci) dfci->compute();
      else assert(false);
      result = fci ? fci->energy() : dfci->energy();
#ifndef DISABLE_SERIALIZATION
    } else if (method == "continue") {
      IArchive archive(itree->get<std::string>("archive"));
      Method* ptr;
      archive >> ptr;
      auto fci = std::shared_ptr<Method>(ptr);
      fci->compute();
      result = std::dynamic_pointer_cast<FCI>(fci)->energy();
#endif
    } else {
      throw std::logic_error("unknown method in fci test");
    }
  }
  std::cout.rdbuf(backup_stream);
  return result;
}

std::vector<double> reference_fci_energy() {
  std::vector<double> out(2);
  out[0] = -98.56280393;
  out[1] = -98.36567235;
  return out;
}
std::vector<double> reference_fci_energy2() {
  std::vector<double> out(2);
  out[0] = -3.30315793;
  out[1] = -2.78407879;
  return out;
}


BOOST_AUTO_TEST_SUITE(TEST_FCI)

BOOST_AUTO_TEST_CASE(KNOWLES_HANDY) {
    BOOST_CHECK(compare(fci_energy("hf_sto3g_fci_kh"), reference_fci_energy()));
#ifndef DISABLE_SERIALIZATION
//  BOOST_CHECK(compare(fci_energy("hf_sto3g_fci_restart"), reference_fci_energy()));
#endif
    BOOST_CHECK(compare(fci_energy("hhe_svp_fci_kh_trip"), reference_fci_energy2()));
}

BOOST_AUTO_TEST_CASE(HARRISON_ZARRABIAN) {
    BOOST_CHECK(compare(fci_energy("hf_sto3g_fci_hz"), reference_fci_energy()));
    BOOST_CHECK(compare(fci_energy("hhe_svp_fci_hz_trip"), reference_fci_energy2()));
}

#ifdef HAVE_MPI_H
BOOST_AUTO_TEST_CASE(DIST_FCI) {
    BOOST_CHECK(compare(fci_energy("hf_sto3g_fci_dist"), reference_fci_energy()));
}
#endif

BOOST_AUTO_TEST_SUITE_END()

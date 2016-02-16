//
// BAGEL - Parallel electron correlation program.
// Filename: test_zfs.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/multi/zcasscf/zcasbfgs.h>
#include <src/prop/pseudospin/pseudospin.h>

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

    } else if (method == "zcasscf") {
      std::string algorithm = itree->get<std::string>("algorithm", "");
      if (algorithm == "bfgs") {
        auto zcas = std::make_shared<ZCASBFGS>(itree, geom, ref);
        zcas->compute();
        ref = zcas->conv_to_ref();
        std::shared_ptr<const PTree> aniso_data = itree->get_child_optional("aniso");
        if (aniso_data) {
          const int nspin = aniso_data->get<int>("nspin", 3);
          auto ps = std::make_shared<Pseudospin>(nspin, aniso_data);
          ps->compute(*zcas->fci());
          out.clear();
          out.push_back(ps->gval(0));
          out.push_back(ps->gval(1));
          out.push_back(ps->gval(2));
          out.push_back(ps->Dval());
          out.push_back(ps->Eval());
        }
      }
    }
  }
  std::cout.rdbuf(backup_stream);
  return out;
}

std::vector<double> reference_pseudospin_parameters() {
  std::vector<double> out(5);
  out[0] =  4.00860253;
  out[1] =  4.00860237;
  out[2] =  4.00442919;
  out[3] =  0.00000164;
  out[4] =  0.00000000;
  return out;
}

BOOST_AUTO_TEST_SUITE(TEST_PSEUDOSPIN)

BOOST_AUTO_TEST_CASE(ZFSPLITTING) {
  BOOST_CHECK(compare(zfs_param("o2_sto3g_zcasscf_pseudospin"), reference_pseudospin_parameters()));
}

BOOST_AUTO_TEST_SUITE_END()

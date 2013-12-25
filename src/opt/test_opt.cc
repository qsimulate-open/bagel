//
// BAGEL - Parallel electron correlation program.
// Filename: test_opt.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/opt/optimize.h>
#include <src/scf/scf.h>
#include <src/wfn/reference.h>

std::vector<double> run_opt(std::string filename) {

  std::string outputname = filename + ".testout";
  std::string inputname = "../../test/" + filename + ".json";
  auto ofs = std::make_shared<std::ofstream>(outputname, std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  auto idata = std::make_shared<const PTree>(inputname);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<const Geometry> geom;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<const Geometry>(itree);
    } else {
      auto opt = std::make_shared<Optimize>(itree, geom);
      opt->compute();

      std::cout.rdbuf(backup_stream);
      std::shared_ptr<const Matrix> out = opt->geometry()->xyz();
      return std::vector<double>(out->data(), out->data()+out->size());
    }
  }
  assert(false);
  return std::vector<double>();
}
std::vector<double> reference_scf_opt() {
  std::vector<double> out(6);
  out[2] = 1.749334;
  out[5] = 0.047492;
  return out;
}
std::vector<double> reference_scf_opt_cart() {
  std::vector<double> out(6);
  out[2] = 1.719396;
  out[5] = 0.016940;
  return out;
}
std::vector<double> reference_scf_opt_mix() {
  std::vector<double> out(6);
  out[2] = 1.755386;
  out[5] =-0.006972;
  return out;
}
std::vector<double> reference_uhf_opt() {
  std::vector<double> out(6);
  out[2] = 1.800736;
  out[5] =-0.005890;
  return out;
}
std::vector<double> reference_rohf_opt() {
  std::vector<double> out(6);
  out[2] = 2.014163;
  out[5] =-0.084977;
  return out;
}
std::vector<double> reference_ks_opt() {
  std::vector<double> out(6);
  out[2] = 1.749755;
  out[5] = 0.002208;
  return out;
}
std::vector<double> reference_cas_act_opt() {
  std::vector<double> out(6);
  out[2] = 1.734489;
  out[5] =-0.003832;
  return out;
}
std::vector<double> reference_mp2_opt() {
  std::vector<double> out(6);
  out[2] = 1.932841;
  out[5] = 0.195935;
  return out;
}
std::vector<double> reference_mp2_aux_opt() {
  std::vector<double> out(6);
  out[2] = 1.932841;
  out[5] = 0.195933;
  return out;
}
std::vector<double> reference_dcf_opt() {
  std::vector<double> out(6);
  out[2] = 0.187050;
  out[5] =-1.532181;
  return out;
}
std::vector<double> reference_hcl_opt() {
  std::vector<double> out(6);
  out[2] = 2.719961;
  out[5] = 0.317095;
  return out;
}

BOOST_AUTO_TEST_SUITE(TEST_OPT)

BOOST_AUTO_TEST_CASE(DF_HF_Opt) {
    BOOST_CHECK(compare(run_opt("hf_svp_dfhf_opt"),       reference_scf_opt(),      1.0e-4));
    BOOST_CHECK(compare(run_opt("hf_svp_dfhf_opt_cart"),  reference_scf_opt_cart(), 1.0e-4));
    BOOST_CHECK(compare(run_opt("hf_mix_dfhf_opt"),       reference_scf_opt_mix(),  1.0e-4));
    BOOST_CHECK(compare(run_opt("oh_svp_uhf_opt"),        reference_uhf_opt(),      1.0e-4));
    BOOST_CHECK(compare(run_opt("hc_svp_rohf_opt"),       reference_rohf_opt(),     1.0e-4));
    BOOST_CHECK(compare(run_opt("hcl_svp_dfhf_opt"),      reference_hcl_opt(),      1.0e-4));
    BOOST_CHECK(compare(run_opt("hf_svp_coulomb_opt"),    reference_dcf_opt(),      1.0e-4));
}
#ifdef HAVE_XC_H
BOOST_AUTO_TEST_CASE(DF_KS_Opt) {
    BOOST_CHECK(compare(run_opt("hf_svp_b3lyp_opt"),      reference_ks_opt(),       1.0e-4));
}
#endif
BOOST_AUTO_TEST_CASE(MP2_Opt) {
    BOOST_CHECK(compare<std::vector<double>>(run_opt("hf_svp_mp2_opt"),        reference_mp2_opt(),      1.0e-4));
    BOOST_CHECK(compare<std::vector<double>>(run_opt("hf_svp_mp2_aux_opt"),    reference_mp2_aux_opt(),  1.0e-4));
}
BOOST_AUTO_TEST_CASE(CASSCF_Opt) {
    BOOST_CHECK(compare<std::vector<double>>(run_opt("hf_svp_cas_act_opt"),    reference_cas_act_opt(),      1.0e-4));
}

BOOST_AUTO_TEST_SUITE_END()

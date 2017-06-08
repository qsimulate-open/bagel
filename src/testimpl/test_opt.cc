//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_opt.cc
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

#include <src/opt/optimize.h>
#include <src/wfn/reference.h>

std::vector<double> run_opt(std::string filename) {

  std::string outputname = filename + ".testout";
  std::string inputname = location__ + filename + ".json";
  auto ofs = std::make_shared<std::ofstream>(outputname, std::ios::trunc);
  std::streambuf* backup_stream = std::cout.rdbuf(ofs->rdbuf());

  auto idata = std::make_shared<const PTree>(inputname);
  auto keys = idata->get_child("bagel");
  std::shared_ptr<const Geometry> geom;
  std::shared_ptr<const Reference> ref;

  std::vector<double> out;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<const Geometry>(itree);
    } else {
      auto opt = std::make_shared<Optimize>(itree, geom, ref);
      opt->compute();

      std::shared_ptr<const Matrix> tmp = opt->geometry()->xyz();
      out = std::vector<double>(tmp->data(), tmp->data()+tmp->size());
    }
  }
  assert(!out.empty());
  std::cout.rdbuf(backup_stream);
  return out;
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
std::vector<double> reference_cas_opt() {
  std::vector<double> out(6);
  out[2] = 1.734489;
  out[5] =-0.003832;
  return out;
}
std::vector<double> reference_sacas_opt() {
  std::vector<double> out(6);
  out[2] = 1.702348;
  out[5] =-0.000261;
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
std::vector<double> reference_ch2_opt() {
  std::vector<double> out(9);
  out[0] =-7.922434;
  out[1] = 5.294612;
  out[2] =-0.676026;
  out[3] =-8.364093;
  out[4] = 3.365386;
  out[5] =-0.356820;
  out[6] =-7.482149;
  out[7] = 7.223682;
  out[8] =-0.997476;
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
    BOOST_CHECK(compare<std::vector<double>>(run_opt("hf_svp_cas_opt"),    reference_cas_opt(),      1.0e-4));
    BOOST_CHECK(compare<std::vector<double>>(run_opt("hf_svp_sacas_opt"),  reference_sacas_opt(),    1.0e-4));
    BOOST_CHECK(compare<std::vector<double>>(run_opt("ch2_sto3g_meci_opt"),reference_ch2_opt(),      1.0e-4));
}

BOOST_AUTO_TEST_SUITE_END()

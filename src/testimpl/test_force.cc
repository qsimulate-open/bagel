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

#include <src/grad/force.h>
#include <src/wfn/reference.h>

std::vector<double> run_force(std::string filename) {

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
    } else if (method == "force") {
      auto force = std::make_shared<Force>(itree, geom, ref);
      std::shared_ptr<const GradFile> grad = force->compute();
      out = std::vector<double>(grad->data(), grad->data()+grad->size());
    } else if (method == "nacme") {
      // sign of nacme depends on the phase of wfn, so compare absolute values
      auto force = std::make_shared<Force>(itree, geom, ref);
      std::shared_ptr<const GradFile> grad = force->compute();
      out = std::vector<double>(grad->data(), grad->data()+grad->size());
      for (auto& i : out) i = fabs(i);
    } else {
      throw std::logic_error("Not yet implemented (run_force)");
    }
  }
  assert(!out.empty());
  std::cout.rdbuf(backup_stream);
  return out;
}
std::vector<double> reference_scf_finite_mix() {
  std::vector<double> out(6);
  out[2] = -0.416074;
  out[5] =  0.416074;
  return out;
}
std::vector<double> reference_svp_mp2_aux_finite() {
  std::vector<double> out(6);
  out[2] = -0.100207;
  out[5] =  0.100207;
  return out;
}
std::vector<double> reference_xms_finite() {
  std::vector<double> out(6);
  out[2] =  0.0336409542;
  out[5] = -0.0336409542;
  return out;
}
std::vector<double> reference_xms_nacme() {
  std::vector<double> out(6);
  out[2] =  0.1722559726;
  out[5] =  0.0786702782;
  return out;
}
std::vector<double> reference_cas_nacme() {
  std::vector<double> out(6);
  out[2] =  0.2075313638;
  out[5] =  0.1243972419;
  return out;
}
std::vector<double> reference_dkh_grad() {
  std::vector<double> out(6);
  out[2] = 0.1590879594;
  out[5] = -0.1590879594;
  return out;
}

BOOST_AUTO_TEST_SUITE(TEST_FORCE)

BOOST_AUTO_TEST_CASE(Deriv_Coup) {
    BOOST_CHECK(compare(run_force("lif_svp_cas_nacme"),      reference_cas_nacme(), 1.0e-5));
#ifdef COMPILE_SMITH
    BOOST_CHECK(compare(run_force("lif_svp_xmscaspt2_nacme"), reference_xms_nacme(), 1.0e-5));
#endif
}

BOOST_AUTO_TEST_CASE(Finite_Grad) {
    BOOST_CHECK(compare(run_force("hf_mix_dfhf_finite"),     reference_scf_finite_mix(), 1.0e-5));
    BOOST_CHECK(compare(run_force("hf_svp_mp2_aux_finite"),  reference_svp_mp2_aux_finite(), 1.0e-5));
#ifdef COMPILE_SMITH
    BOOST_CHECK(compare(run_force("lif_svp_xmscaspt2_finite"), reference_xms_finite(), 1.0e-5));
#endif
}

BOOST_AUTO_TEST_CASE(Hcore_Grad) {
    BOOST_CHECK(compare(run_force("hf_svp_dfhf_dkh_grad"),   reference_dkh_grad(), 1.0e-5));
}

BOOST_AUTO_TEST_SUITE_END()

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_smith.cc
// Copyright (C) 2016 Toru Shiozaki
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

#include <src/wfn/reference.h>

std::vector<double> reference_noshift() {
  std::vector<double> out(6);
  out[2] =  0.0043773889;
  out[5] = -0.0043773889;
  return out;
}
std::vector<double> reference_shift() {
  std::vector<double> out(6);
  out[2] =  0.0042621205;
  out[5] = -0.0042621205;
  return out;
}
std::vector<double> reference_ms() {
  std::vector<double> out(6);
  out[2] =  0.0396123988;
  out[5] = -0.0396123988;
  return out;
}
std::vector<double> reference_xms() {
  std::vector<double> out(6);
  out[2] =  0.0336409542;
  out[5] = -0.0336409542;
  return out;
}
std::vector<double> reference_xms_imag() {
  std::vector<double> out(6);
  out[2] =  0.0340340025;
  out[5] = -0.0340340025;
  return out;
}

#ifdef COMPILE_SMITH
BOOST_AUTO_TEST_SUITE(TEST_SMITH)

BOOST_AUTO_TEST_CASE(CASPT2_Opt) {
    BOOST_CHECK(compare(run_force("li2_svp_caspt2_grad"),    reference_noshift(),  1.0e-5));
    BOOST_CHECK(compare(run_force("li2_svp_caspt2_shift"),   reference_shift(),  1.0e-5));
    BOOST_CHECK(compare(run_force("lif_svp_mscaspt2_grad"),  reference_ms(),  1.0e-5));
    BOOST_CHECK(compare(run_force("lif_svp_xmscaspt2_grad"), reference_xms(), 1.0e-5));
    BOOST_CHECK(compare(run_force("lif_svp_xmscaspt2_grad_imag"), reference_xms_imag(), 1.0e-5));
}

BOOST_AUTO_TEST_SUITE_END()
#endif

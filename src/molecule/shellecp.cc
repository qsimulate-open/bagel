//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: shell_ecp.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <sstream>
#include <src/molecule/shellecp.h>

using namespace std;
using namespace bagel;

Shell_ECP::Shell_ECP()
 : Shell_base(false, {{0.0, 0.0, 0.0}}, 0),
   ecp_exponents_(0.0), ecp_coefficients_(0.0), ecp_r_power_(0) {}

Shell_ECP::Shell_ECP(const array<double,3>& _position, const int _ang,
                     const vector<double>& _ecp_expo, const vector<double>& _ecp_coef,
                     const vector<int>& _ecp_r)
 : Shell_base(false, _position, _ang),
   ecp_exponents_(_ecp_expo), ecp_coefficients_(_ecp_coef), ecp_r_power_(_ecp_r) {}

string Shell_ECP::show() const {
  stringstream ss;
  ss << "position: ";
  ss << position_[0] << " " << position_[1] << " "  << position_[2] << endl;
  ss << "ecp_angular: "  << angular_number_ << endl;
  ss << "ecp_exponents: ";
  for (int i = 0; i != ecp_exponents_.size(); ++i) {
    ss << " " << ecp_exponents_[i];
  }
  ss << endl;
  ss << "ecp_coefficients: ";
  for (int i = 0; i != ecp_coefficients_.size(); ++i) {
      ss << " " << ecp_coefficients_[i];
  }
  ss << endl;
  ss << "ecp_r_power: ";
  for (int i = 0; i != ecp_r_power_.size(); ++i) {
      ss << " " << ecp_r_power_[i];
  }
  ss << endl;

  return ss.str();
}


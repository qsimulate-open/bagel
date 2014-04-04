//
// BAGEL - Parallel electron correlation program.
// Filename: shell_ECP.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/molecule/shell_ECP.h>

using namespace std;
using namespace bagel;

Shell_ECP::Shell_ECP(const bool sph, const array<double,3>& _position, int _ang, const vector<double>& _expo,
                     const vector<vector<double>>& _contr,  const vector<pair<int, int>>& _range,
                     const std::vector<double>& _ecp_expo, const std::vector<double>& _ecp_coef,
                     const std::vector<int>& _ecp_r)
 : Shell_base(sph, _position, _ang, _expo, _contr, _range),
   ecp_exponents_(_ecp_expo), ecp_coefficients_(_ecp_coef), ecp_r_power_(_ecp_r) {}


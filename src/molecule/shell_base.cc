//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: shell_base.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <sstream>
#include <src/molecule/shell_base.h>

using namespace std;
using namespace bagel;

Shell_base::Shell_base(const bool sph, const array<double,3>& _position, int _ang)
 : spherical_(sph), position_(_position), angular_number_(_ang) {}

Shell_base::Shell_base(const bool sph) : spherical_(sph), position_{{0.0,0.0,0.0}}, angular_number_(0) {}

std::string Shell_base::show() const {
  stringstream ss;
  ss << "position: ";
  ss << position_[0] << " " << position_[1] << " "  << position_[2] << endl;

  return ss.str();
}

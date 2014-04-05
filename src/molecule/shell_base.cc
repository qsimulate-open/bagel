//
// BAGEL - Parallel electron correlation program.
// Filename: shell_base.cc
// Copyright (C) 2009 Toru Shiozaki
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

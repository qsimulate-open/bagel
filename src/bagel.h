//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: bagel.h
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


#ifndef __SRC_BAGEL__H
#define __SRC_BAGEL__H

#include <string>

namespace bagel {

// run BAGEL from an input file
extern void run_bagel_from_input(const std::string& filename);
// run BAGEL from stringify'ed json
extern void run_bagel_from_json(const std::string& json);

}

#endif


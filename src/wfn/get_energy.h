//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: get_energy.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_WFN_GET_ENERGY_H
#define __SRC_WFN_GET_ENERGY_H

#include <src/wfn/method.h>
#include <src/wfn/reference.h>
#include <src/util/input/input.h>

namespace bagel {
  extern std::tuple<double, std::shared_ptr<const Reference>> get_energy(std::string title, std::shared_ptr<const PTree> itree,
            std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref, const int target = 0);
}

#endif

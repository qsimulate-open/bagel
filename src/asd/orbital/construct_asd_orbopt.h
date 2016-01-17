//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd/orbital/construct_asd_orbopt.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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


#ifndef __SRC_ASD_CONSTRUCT_ORBOPT_H
#define __SRC_ASD_CONSTRUCT_ORBOPT_H

#include <src/asd/orbital/asd_orbopt.h>
#include <src/asd/dimer/dimer.h>

namespace bagel {
  extern std::shared_ptr<ASD_OrbOpt> construct_ASD_OrbOpt(std::shared_ptr<const PTree> itree, std::shared_ptr<Dimer> dimer);
}

#endif

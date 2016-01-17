//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: construct_asd.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_ASD_CONSTRUCT_ASD_H
#define __SRC_ASD_CONSTRUCT_ASD_H

#include <src/asd/asd_base.h>
#include <src/asd/dimer/dimer.h>

namespace bagel {
  extern std::shared_ptr<ASD_base> construct_ASD(std::shared_ptr<const PTree> itree, std::shared_ptr<Dimer> dimer, bool rdm = false);
}

#endif

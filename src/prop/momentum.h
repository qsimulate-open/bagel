//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: momentum.h
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


#ifndef __SRC_PROP_MOMENTUM_H
#define __SRC_PROP_MOMENTUM_H

#include <src/wfn/geometry.h>

namespace bagel {

class Momentum {
  protected:
    std::shared_ptr<const Geometry> geom_;

  public:
    Momentum(std::shared_ptr<const Geometry>);

    std::array<std::shared_ptr<Matrix>,3> compute() const;
};

}

#endif

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: momentum_point.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_PROP_MOMENTUM_POINT_H
#define __SRC_PROP_MOMENTUM_POINT_H

#include <src/wfn/geometry.h>

namespace bagel {

class Momentum_Point {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::array<double,3> location_;

  public:
    Momentum_Point(std::shared_ptr<const Geometry> g, std::array<double,3> loc);

    std::array<std::shared_ptr<ZMatrix>,3> compute() const;
};

}

#endif

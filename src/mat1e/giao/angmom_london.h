//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: angmom_london.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __SRC_MAT1E_GIAO_ANGMOM_LONDON_H
#define __SRC_MAT1E_GIAO_ANGMOM_LONDON_H

#include <src/wfn/geometry.h>

namespace bagel {

class AngMom_London {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::array<double,3> mcoord_;

  public:
    AngMom_London(std::shared_ptr<const Geometry> geom, std::array<double,3> mcoord);

    std::array<std::shared_ptr<ZMatrix>,3> compute() const;
};

}

#endif

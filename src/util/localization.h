//
// BAGEL - Parallel electron correlation program.
// Filename: localization.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __BAGEL_UTIL_LOCALIZE_H
#define __BAGEL_UTIL_LOCALIZE_H

#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <src/scf/geometry.h>
#include <src/util/matrix.h>

namespace bagel {

class OrbitalLocalization {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<Matrix> density_;

  public:
    OrbitalLocalization(std::shared_ptr<const Geometry> geom, std::shared_ptr<Matrix> density) : geom_(geom), density_(density) {};

    virtual std::shared_ptr<Matrix> localize() = 0;
};

class RegionLocalization : public OrbitalLocalization {
  protected:
    std::vector<std::pair<int, int> > bounds_;

  public:
    RegionLocalization(std::shared_ptr<const Geometry> geom, std::shared_ptr<Matrix> density, std::vector<std::pair<int, int> > atom_bounds);

    std::shared_ptr<Matrix> localize() override;
  
};

}

#endif

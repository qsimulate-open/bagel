//
// BAGEL - Parallel electron correlation program.
// Filename: spinint.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

// Produces matrix elements of S_x, S_y, and S_z in atomic orbital basis

#ifndef __SRC_MAT1E_REL_SPININT_H
#define __SRC_MAT1E_REL_SPININT_H

#include <src/wfn/geometry.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

class RelSpinInt {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    void compute_();

    std::array<std::shared_ptr<ZMatrix>,3> data_;

  public:
    RelSpinInt(const std::shared_ptr<const Geometry> geom) : geom_(geom) {
      assert(!geom->magnetism()); // GIAO-RMB version has not been implemented
      compute_();
    }

    std::shared_ptr<ZMatrix>operator()(const int i) const { assert(i >= 0 && i < 3); return data_[i]; }
};

}

#endif

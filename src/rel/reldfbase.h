//
// BAGEL - Parallel electron correlation program.
// Filename: dfhalfbase.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#ifndef __SRC_REL_RELDFBASE_H
#define __SRC_REL_RELDFBASE_H

#include <src/rel/alpha.h>
#include <src/math/zmatrix.h>
#include <src/rel/spinorinfo.h>

namespace bagel {

class RelDFBase {
  protected:
    // l,x,y,z
    std::pair<int, int> cartesian_;
    // X,Y, and coefficients
    std::vector<std::shared_ptr<const SpinorInfo>> basis_;

  public:
    RelDFBase(std::pair<int, int> cartesian) : cartesian_(cartesian) {
    }

    RelDFBase(const RelDFBase& o) : cartesian_(o.cartesian_) {
    }

    std::pair<int, int> cartesian() const { return cartesian_; }
    const std::vector<std::shared_ptr<const SpinorInfo>>& basis() const { return basis_; }

};

}

#endif

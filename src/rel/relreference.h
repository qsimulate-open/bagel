//
// BAGEL - Parallel electron correlation program.
// Filename: relreference.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_REL_RELREFERENCE_H
#define __SRC_REL_RELREFERENCE_H

#include <src/wfn/geometry.h>
#include <src/util/zmatrix.h>

namespace bagel {

class RelReference {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const ZMatrix> coeff_; 

  public:
    RelReference(std::shared_ptr<const Geometry> g, std::shared_ptr<const ZMatrix> c)
     : geom_(g), coeff_(c) {
    }

    const std::shared_ptr<const Geometry> geom() const { return geom_; }
    const std::shared_ptr<const ZMatrix> coeff() const { return coeff_; }

};

}

#endif

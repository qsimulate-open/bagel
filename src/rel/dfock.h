//
// BAGEL - Parallel electron correlation program.
// Filename: dfock.h
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


#ifndef __SRC_REL_DFOCK_H
#define __SRC_REL_DFOCK_H

#include <memory>
#include <string>
#include <map>
#include <src/wfn/reference.h>
#include <src/scf/geometry.h>
#include <src/util/zmatrix.h>

namespace bagel {

class DFock : public ZMatrix {
  protected:
    std::shared_ptr<const Geometry> geom_;
    void two_electron_part(const std::shared_ptr<const ZMatrix> ocoeff, const bool rhf, const double scale_ex);

  public:
    DFock(const std::shared_ptr<const Geometry> a, 
          const std::shared_ptr<const ZMatrix> ocoeff, const bool rhf = false, const double scale_ex = 1.0)
     : ZMatrix(a->nbasis(), a->nbasis()), geom_(a) {
       two_electron_part(ocoeff, rhf, scale_ex);
    }

//    std::shared_ptr<Reference> conv_to_ref() const override;

};

}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: sofock.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __BAGEL_SRC_SCF_SOFOCK_H
#define __BAGEL_SRC_SCF_SOFOCK_H

#include <src/df/df.h>
#include <src/wfn/geometry.h>
#include <src/integral/libint/libint.h>
#include <src/math/zmatrix.h>

namespace bagel {

class SOFock : public ZMatrix {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const ZMatrix> previous_;
    std::shared_ptr<const ZMatrix> socoeff_;
    void form_sofock();

  public:
    SOFock(const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const ZMatrix> previous, const std::shared_ptr<const ZMatrix> socoeff);

};

}
#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: gkineticbatch.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_INTEGRAL_OS_GKINETICBATCH_H
#define __SRC_INTEGRAL_OS_GKINETICBATCH_H

#include <src/integral/os/osintegral.h>

namespace bagel {

// computes derivative integrals of kinetic operator.
class GKineticBatch : public OSInt {
  protected:
    int nblocks() const override { return 6; }
    int nrank() const override { return 3; }

  public:
    GKineticBatch(const std::array<std::shared_ptr<const Shell>,2>& o) : OSInt(o) { common_init(); }

    void compute();

};

}

#endif

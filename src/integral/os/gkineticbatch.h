//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gkineticbatch.h
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

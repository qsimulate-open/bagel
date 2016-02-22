//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: overlapbatch.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_INTEGRAL_OS_OVERLAPBATCH_H
#define __SRC_INTEGRAL_OS_OVERLAPBATCH_H

#include <src/integral/os/osintegral.h>

namespace bagel {

class OverlapBatch : public OSInt {
  protected:
    void perform_VRR(double*) override;

    int nblocks() const override { return 1; }
    int nrank() const override { return 0; }

  public:
    OverlapBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, std::shared_ptr<StackMem> stack = nullptr)
    : OSInt(basis, stack) { common_init(); }

    constexpr static int Nblocks() { return 1; }
    void compute() override;
};

}

#endif

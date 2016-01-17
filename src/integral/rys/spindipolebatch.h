//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spindipole.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_INTEGRAL_RYS_SPINDIPOLEBATCH_H
#define __SRC_INTEGRAL_RYS_SPINDIPOLEBATCH_H

#include <src/integral/rys/coulombbatch_base.h>

namespace bagel {

class SpinDipoleBatch : public CoulombBatch_Base<double> {

  protected:
    std::shared_ptr<const Atom> target_;

    void root_weight(const int ps) override;

  public:

    SpinDipoleBatch(const std::array<std::shared_ptr<const Shell>,2>& _info, std::shared_ptr<const Atom> target,
                    std::shared_ptr<StackMem> stack = nullptr);

    void compute() override;

    double* data(const int i) override { return data_ + i*size_block_; }

    constexpr static int Nblocks() { return 6; }

};

}

#endif

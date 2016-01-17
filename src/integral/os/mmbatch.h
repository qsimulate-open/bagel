//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mmbatch.h
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


#ifndef __SRC_INTEGRAL_OS_MMBATCH_H
#define __SRC_INTEGRAL_OS_MMBATCH_H

#include <src/integral/os/osintegral.h>
#include <src/molecule/molecule.h>

namespace bagel {

class MMBatch : public OSInt {
  protected:
    const std::array<double,3> center_;
    void perform_VRR(double*) override;

    int nblocks() const override { return lmax_*(lmax_*lmax_+6*lmax_+11)/6; } // \sum_1^n (j+1)(j+2)/2
    int nrank() const override { return 0; }

    const int lmax_;

  public:
    MMBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, std::shared_ptr<const Molecule> c, const int lmax)
     : OSInt(basis), center_(c->charge_center()), lmax_(lmax) { common_init(); }

    void compute() override;

    int num_blocks() const { return nblocks(); }

};

class DipoleBatch : public MMBatch {
  public:
    DipoleBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, std::shared_ptr<const Molecule> c) : MMBatch(basis, c, 1)
      { assert(nblocks() == Nblocks()); }
    constexpr static int Nblocks() { return 3; }
};

}

#endif

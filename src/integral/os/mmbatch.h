//
// BAGEL - Parallel electron correlation program.
// Filename: mmbatch.h
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

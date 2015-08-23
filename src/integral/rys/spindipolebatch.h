//
// BAGEL - Parallel electron correlation program.
// Filename: spindipole.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_INTEGRAL_RYS_SPINDIPOLEBATCH_H
#define __SRC_INTEGRAL_RYS_SPINDIPOLEBATCH_H

#include <src/integral/rys/coulombbatch_base.h>

namespace bagel {

class SpinDipoleBatch : public CoulombBatch_Base<double> {

  protected:
    std::shared_ptr<const Atom> target_;

    void root_weight(const int ps) override;
    virtual double scale_root(const double root, const double p, const double zeta) { return root; }
    virtual double scale_weight(const double weight) { return weight; }

  public:

    SpinDipoleBatch(const std::array<std::shared_ptr<const Shell>,2>& _info, std::shared_ptr<const Atom> target,
                    std::shared_ptr<StackMem> stack = nullptr)
      : CoulombBatch_Base<double>(_info, std::make_shared<Molecule>(std::vector<std::shared_ptr<const Atom>>{target}, std::vector<std::shared_ptr<const Atom>>{}),
                                  0, 2, stack), target_(target) {
    }

    void compute() override;

    constexpr static int Nblocks() { return 6; }

};

}

#endif

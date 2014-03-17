//
// BAGEL - Parallel electron correlation program.
// Filename: coulombbatch_energy.h
// Copyright (C) 2009 Toru Shiozaki
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

#ifndef __SRC_INTEGRAL_RYS_COULOMBBATCH_ENERGY_H
#define __SRC_INTEGRAL_RYS_COULOMBBATCH_ENERGY_H

#include <src/integral/rys/coulombbatch_base.h>

namespace bagel {

class CoulombBatch_energy : public CoulombBatch_Base<double> {

  protected:

    void root_weight(const int ps) override;
    virtual double scale_root(const double root, const double p, const double zeta) { return root; }
    virtual double scale_weight(const double weight, const double coef) { return weight; }

  public:

    CoulombBatch_energy(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, std::shared_ptr<StackMem> stack = nullptr)
      :  CoulombBatch_Base<double>(_info, mol, 0, stack, 0, 0.0) {};

    CoulombBatch_energy(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, const int L, const double A = 0.0)
      :  CoulombBatch_Base<double>(_info, mol, 0, nullptr, L, A) {};

     ~CoulombBatch_energy() {};

    void compute() override;

    constexpr static int Nblocks() { return 1; }

};

}

#endif

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

#include <src/integral/sortlist.h>
#include <src/integral/carsphlist.h>
#include <src/integral/hrrlist.h>
#include <src/integral/rys/coulombbatch_energy.h>
#include <src/integral/rys/coulombbatch_base.h>

namespace bagel {

template <typename DataType>
class CoulombBatch_Energy : public CoulombBatch_Base<DataType> {

  protected:

    virtual double scale_root(const double root, const double p, const double zeta) { return root; }
    virtual double scale_weight(const double weight, const double coef) { return weight; }

  public:

    CoulombBatch_Energy(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>())
      :  CoulombBatch_Base<DataType>(_info, mol, 0, stack, 0, 0.0) {};

    CoulombBatch_Energy(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, const int L, const double A = 0.0)
      :  CoulombBatch_Base<DataType>(_info, mol, 0, std::shared_ptr<StackMem>(), L, A) {};

     ~CoulombBatch_Energy() {};

    void compute() override;

    constexpr static int Nblocks() { return 1; }

  protected:

    using CoulombBatch_Base<DataType>::natom_;
    using CoulombBatch_Base<DataType>::L_;
    using CoulombBatch_Base<DataType>::A_;
    using CoulombBatch_Base<DataType>::mol_;

    using RysIntegral<DataType>::bkup_;
    using RysIntegral<DataType>::spherical1_;
    using RysIntegral<DataType>::basisinfo_;
    using RysIntegral<DataType>::prim0size_;
    using RysIntegral<DataType>::prim1size_;
    using RysIntegral<DataType>::contsize_;
    using RysIntegral<DataType>::cont0size_;
    using RysIntegral<DataType>::cont1size_;
    using RysIntegral<DataType>::swap01_;
    using RysIntegral<DataType>::screening_;
    using RysIntegral<DataType>::screening_size_;
    using RysIntegral<DataType>::size_alloc_;
    using RysIntegral<DataType>::size_final_;
    using RysIntegral<DataType>::stack_;
    using RysIntegral<DataType>::xp_;
    using RysIntegral<DataType>::P_;
    using RysIntegral<DataType>::AB_;
    using RysIntegral<DataType>::coeff_;
    using RysIntegral<DataType>::rank_;
    using RysIntegral<DataType>::roots_;
    using RysIntegral<DataType>::weights_;
    using RysIntegral<DataType>::asize_;
    using RysIntegral<DataType>::amapping_;
    using RysIntegral<DataType>::amax_;
    using RysIntegral<DataType>::amin_;
    using RysIntegral<DataType>::amax1_;
    using RysIntegral<DataType>::data_;

};

using CoulombBatch_energy = CoulombBatch_Energy<double>;

}

#define COULOMBBATCH_ENERGY_HEADERS
#include <src/integral/rys/coulombbatch_energy_impl.hpp>
#undef COULOMBBATCH_ENERGY_HEADERS

#endif

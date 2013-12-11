//
// BAGEL - Parallel electron correlation program.
// Filename: coulombbatch_base.h
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


#ifndef __SRC_RYSINT_COULOMBBATCH_BASE_H
#define __SRC_RYSINT_COULOMBBATCH_BASE_H

#include <src/integral/rys/rysintegral.h>
#include <src/integral/rys/coulombbatch_base.h>
#include <src/molecule/molecule.h>
#include <src/util/constants.h>
#include <src/integral/rys/inline.h>
#include <src/integral/rys/erirootlist.h>


namespace bagel {

template <typename DataType>
class CoulombBatch_Base : public RysIntegral<DataType> {

  protected:

    std::shared_ptr<const Molecule> mol_;
    int natom_;

    /// for periodic calculations (UNCHECKED!!)
    const int L_;
    const double A_;

    void root_weight(const int ps) override;
    void compute_ssss(const double) override;
    void allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) override;

  public:
    CoulombBatch_Base(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, const int deriv,
                  std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>(),
                  const int L = 0, const double A = 0.0);
    ~CoulombBatch_Base() {}

    std::shared_ptr<const Molecule> mol() const { return mol_; }

  protected:
    using RysIntegral<DataType>::tenno_;
    using RysIntegral<DataType>::breit_;
    using RysIntegral<DataType>::basisinfo_;
    using RysIntegral<DataType>::primsize_;
    using RysIntegral<DataType>::contsize_;
    using RysIntegral<DataType>::screening_;
    using RysIntegral<DataType>::screening_size_;
    using RysIntegral<DataType>::size_block_;
    using RysIntegral<DataType>::size_alloc_;
    using RysIntegral<DataType>::size_final_;
    using RysIntegral<DataType>::stack_;
    using RysIntegral<DataType>::stack_save_;
    using RysIntegral<DataType>::stack_save2_;
    using RysIntegral<DataType>::xp_;
    using RysIntegral<DataType>::P_;
    using RysIntegral<DataType>::AB_;
    using RysIntegral<DataType>::T_;
    using RysIntegral<DataType>::coeff_;
    using RysIntegral<DataType>::rank_;
    using RysIntegral<DataType>::roots_;
    using RysIntegral<DataType>::weights_;
    using RysIntegral<DataType>::deriv_rank_;
    using RysIntegral<DataType>::asize_;
    using RysIntegral<DataType>::csize_;
    using RysIntegral<DataType>::amax_;
    using RysIntegral<DataType>::cmax_;
    using RysIntegral<DataType>::data_;
    using RysIntegral<DataType>::data2_;

};

using CoulombBatch_base = CoulombBatch_Base<double>;

}

#define COULOMBBATCH_BASE_HEADERS
#include <src/integral/rys/coulombbatch_base_impl.hpp>
#undef COULOMBBATCH_BASE_HEADERS

#endif

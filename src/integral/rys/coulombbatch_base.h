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
#include <src/molecule/molecule.h>
#include <src/util/constants.h>
#include <src/integral/rys/inline.h>
#include <src/integral/rys/erirootlist.h>


namespace bagel {

template <typename DataType, Int_t IntType = Int_t::Standard>
class CoulombBatch_Base : public RysIntegral<DataType, IntType> {

  protected:

    std::shared_ptr<const Molecule> mol_;
    int natom_;

    void compute_ssss(const double) override;
    void allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) override;
    virtual DataType get_PQ(const double coord1, const double coord2, const double exp1, const double exp2, const double one12, const int center1, const int dim, const bool swap) {
      return (coord1*exp1 + coord2*exp2) * one12;
    }

  public:
    CoulombBatch_Base(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, const int deriv,
                  std::shared_ptr<StackMem> stack = nullptr);
    ~CoulombBatch_Base() {}

    std::shared_ptr<const Molecule> mol() const { return mol_; }

  protected:
    using RysIntegral<DataType,IntType>::tenno_;
    using RysIntegral<DataType,IntType>::breit_;
    using RysIntegral<DataType,IntType>::basisinfo_;
    using RysIntegral<DataType,IntType>::primsize_;
    using RysIntegral<DataType,IntType>::contsize_;
    using RysIntegral<DataType,IntType>::screening_;
    using RysIntegral<DataType,IntType>::screening_size_;
    using RysIntegral<DataType,IntType>::size_block_;
    using RysIntegral<DataType,IntType>::size_alloc_;
    using RysIntegral<DataType,IntType>::size_final_;
    using RysIntegral<DataType,IntType>::stack_;
    using RysIntegral<DataType,IntType>::stack_save_;
    using RysIntegral<DataType,IntType>::stack_save2_;
    using RysIntegral<DataType,IntType>::xp_;
    using RysIntegral<DataType,IntType>::P_;
    using RysIntegral<DataType,IntType>::AB_;
    using RysIntegral<DataType,IntType>::T_;
    using RysIntegral<DataType,IntType>::coeff_;
    using RysIntegral<DataType,IntType>::rank_;
    using RysIntegral<DataType,IntType>::roots_;
    using RysIntegral<DataType,IntType>::weights_;
    using RysIntegral<DataType,IntType>::deriv_rank_;
    using RysIntegral<DataType,IntType>::asize_;
    using RysIntegral<DataType,IntType>::csize_;
    using RysIntegral<DataType,IntType>::amax_;
    using RysIntegral<DataType,IntType>::cmax_;
    using RysIntegral<DataType,IntType>::data_;
    using RysIntegral<DataType,IntType>::data2_;
    using RysIntegral<DataType,IntType>::swap01_;

};

using CoulombBatch_base = CoulombBatch_Base<double>;

}

#define COULOMBBATCH_BASE_HEADERS
#include <src/integral/rys/coulombbatch_base_impl.hpp>
#undef COULOMBBATCH_BASE_HEADERS

extern template class bagel::CoulombBatch_Base<double>;
//extern template class bagel::CoulombBatch_Base<std::complex<double>>;

#endif

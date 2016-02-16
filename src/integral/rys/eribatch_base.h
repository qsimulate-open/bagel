//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: eribatch_base.h
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


#ifndef __SRC_INTEGRAL_RYS_ERIBATCH_BASE_H
#define __SRC_INTEGRAL_RYS_ERIBATCH_BASE_H

#include <src/integral/rys/inline.h>
#include <src/integral/rys/rysintegral.h>


namespace bagel {

template <typename DataType, Int_t IntType = Int_t::Standard>
class ERIBatch_Base : public RysIntegral<DataType,IntType> {
  protected:
    // a hack for screening of three-center integrals
    static double rnd(const double& a) { return (a > 0.0) ? a : 1.0; };

    virtual void root_weight(const int ps) = 0;

    // this sets T_ (+ U_), P_, Q_, xp_, xq_, coeff_, and screening_size_
    // for ERI evaulation. Other than that, we need to overload this function in a derived class
    void compute_ssss(const double integral_thresh) override;

    virtual DataType get_PQ(const double coord1, const double coord2, const double exp1, const double exp2, const double one12, const int center1, const int dim, const bool swap) {
      return (coord1*exp1 + coord2*exp2) * one12;
    }

    void allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) override;

  public:
    ERIBatch_Base(const std::array<std::shared_ptr<const Shell>,4>& o, const int deriv, const int breit = 0, std::shared_ptr<StackMem> stack = nullptr);

  protected:
    using RysIntegral<DataType,IntType>::basisinfo_;
    using RysIntegral<DataType,IntType>::stack_;
    using RysIntegral<DataType,IntType>::stack_save_;
    using RysIntegral<DataType,IntType>::stack_save2_;
    using RysIntegral<DataType,IntType>::AB_;
    using RysIntegral<DataType,IntType>::CD_;
    using RysIntegral<DataType,IntType>::contsize_;
    using RysIntegral<DataType,IntType>::asize_;
    using RysIntegral<DataType,IntType>::csize_;
    using RysIntegral<DataType,IntType>::primsize_;
    using RysIntegral<DataType,IntType>::prim0size_;
    using RysIntegral<DataType,IntType>::prim1size_;
    using RysIntegral<DataType,IntType>::prim2size_;
    using RysIntegral<DataType,IntType>::prim3size_;
    using RysIntegral<DataType,IntType>::screening_;
    using RysIntegral<DataType,IntType>::screening_size_;
    using RysIntegral<DataType,IntType>::T_;
    using RysIntegral<DataType,IntType>::P_;
    using RysIntegral<DataType,IntType>::Q_;
    using RysIntegral<DataType,IntType>::xp_;
    using RysIntegral<DataType,IntType>::xq_;
    using RysIntegral<DataType,IntType>::size_block_;
    using RysIntegral<DataType,IntType>::size_alloc_;
    using RysIntegral<DataType,IntType>::size_final_;
    using RysIntegral<DataType,IntType>::coeff_;
    using RysIntegral<DataType,IntType>::deriv_rank_;
    using RysIntegral<DataType,IntType>::tenno_;
    using RysIntegral<DataType,IntType>::breit_;
    using RysIntegral<DataType,IntType>::data_;
    using RysIntegral<DataType,IntType>::data2_;
    using RysIntegral<DataType,IntType>::swap01_;
    using RysIntegral<DataType,IntType>::swap23_;
};

using ERIBatch_base = ERIBatch_Base<double, Int_t::Standard>;

}

extern template class bagel::ERIBatch_Base<double,bagel::Int_t::Standard>;
extern template class bagel::ERIBatch_Base<std::complex<double>,bagel::Int_t::London>;

#endif

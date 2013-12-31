//
// BAGEL - Parallel electron correlation program.
// Filename: eribatch_base.h
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


#ifndef __SRC_INTEGRAL_RYS_ERIBATCH_BASE_H
#define __SRC_INTEGRAL_RYS_ERIBATCH_BASE_H

#include <src/integral/rys/inline.h>
#include <src/integral/rys/rysintegral.h>
#include <iomanip>  // Only for testing


namespace bagel {

namespace {
  enum class Int_t { Standard, London };
}

template <typename DataType, Int_t IntType = Int_t::Standard>
class ERIBatch_Base : public RysIntegral<DataType> {

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

    ERIBatch_Base(const std::array<std::shared_ptr<const Shell>,4>& o, const int deriv, const int breit = 0,
              std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>()) : RysIntegral<DataType>(o, stack) {

      breit_ = breit;
      deriv_rank_ = deriv;

      // determines if we want to swap shells
      this->set_swap_info(true);

      // stores AB and CD
      this->set_ab_cd();

      // set primsize_ and contsize_, as well as relevant members
      this->set_prim_contsizes();

      // sets angular info
      int asize_final, csize_final, asize_final_sph, csize_final_sph;
      std::tie(asize_final, csize_final, asize_final_sph, csize_final_sph) = this->set_angular_info();

      // allocate
      allocate_data(asize_final, csize_final, asize_final_sph, csize_final_sph);
      this->allocate_arrays(primsize_);

    }


  protected:
    using RysIntegral<DataType>::basisinfo_;
    using RysIntegral<DataType>::stack_;
    using RysIntegral<DataType>::stack_save_;
    using RysIntegral<DataType>::stack_save2_;
    using RysIntegral<DataType>::AB_;
    using RysIntegral<DataType>::CD_;
    using RysIntegral<DataType>::contsize_;
    using RysIntegral<DataType>::asize_;
    using RysIntegral<DataType>::csize_;
    using RysIntegral<DataType>::primsize_;
    using RysIntegral<DataType>::prim0size_;
    using RysIntegral<DataType>::prim1size_;
    using RysIntegral<DataType>::prim2size_;
    using RysIntegral<DataType>::prim3size_;
    using RysIntegral<DataType>::screening_;
    using RysIntegral<DataType>::screening_size_;
    using RysIntegral<DataType>::T_;
    using RysIntegral<DataType>::P_;
    using RysIntegral<DataType>::Q_;
    using RysIntegral<DataType>::xp_;
    using RysIntegral<DataType>::xq_;
    using RysIntegral<DataType>::size_block_;
    using RysIntegral<DataType>::size_alloc_;
    using RysIntegral<DataType>::size_final_;
    using RysIntegral<DataType>::coeff_;
    using RysIntegral<DataType>::deriv_rank_;
    using RysIntegral<DataType>::tenno_;
    using RysIntegral<DataType>::breit_;
    using RysIntegral<DataType>::data_;
    using RysIntegral<DataType>::data2_;
    using RysIntegral<DataType>::swap01_;
    using RysIntegral<DataType>::swap23_;
};

using ERIBatch_base = ERIBatch_Base<double, Int_t::Standard>;

}

#define ERIBATCH_BASE_HEADERS
#include <src/integral/rys/eribatch_base_impl.hpp>
#undef ERIBATCH_BASE_HEADERS

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: denom.h
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

#ifndef __SRC_SMITH_DENOM_H
#define __SRC_SMITH_DENOM_H

#include <src/wfn/rdm.h>
#include <src/util/math/zmatrix.h>
#include <src/smith/smith_util.h>

namespace bagel {
namespace SMITH {

template<typename DataType>
class Denom {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;

  protected:
    std::shared_ptr<const MatType> fock_;
    const IndexRange active_;
    const IndexRange ortho1_;
    const IndexRange ortho2_;
    const IndexRange ortho3_;
    const IndexRange ortho2t_;
    const double thresh_;

    // work area
    std::shared_ptr<MatType> shalf_x_;
    std::shared_ptr<MatType> shalf_h_;
    std::shared_ptr<MatType> shalf_xx_;
    std::shared_ptr<MatType> shalf_hh_;
    std::shared_ptr<MatType> shalf_xh_;
    std::shared_ptr<MatType> shalf_xhh_;
    std::shared_ptr<MatType> shalf_xxh_;
    // work area
    std::shared_ptr<MatType> work_x_;
    std::shared_ptr<MatType> work_h_;
    std::shared_ptr<MatType> work_xx_;
    std::shared_ptr<MatType> work_hh_;
    std::shared_ptr<MatType> work_xh_;
    std::shared_ptr<MatType> work_xhh_;
    std::shared_ptr<MatType> work_xxh_;

    std::shared_ptr<TATensor<DataType,2>> tashalf_x_;
    std::shared_ptr<TATensor<DataType,2>> tashalf_h_;
    std::shared_ptr<TATensor<DataType,3>> tashalf_xx_;
    std::shared_ptr<TATensor<DataType,3>> tashalf_hh_;
    std::shared_ptr<TATensor<DataType,3>> tashalf_xh_;
    std::shared_ptr<TATensor<DataType,3>> tashalf_xh2_;
    std::shared_ptr<TATensor<DataType,4>> tashalf_xhh_;
    std::shared_ptr<TATensor<DataType,4>> tashalf_xxh_;

    VectorB denom_x_;
    VectorB denom_h_;
    VectorB denom_xx_;
    VectorB denom_hh_;
    VectorB denom_xh_;
    VectorB denom_xhh_;
    VectorB denom_xxh_;

    // init functions
    void init_x_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                       std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_h_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                       std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_xx_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                        std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_hh_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                        std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_xh_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                        std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_xhh_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                         std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    void init_xxh_(const int, const int, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                         std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);

  public:
    Denom(std::shared_ptr<const MatType> fock, const int nstates, const std::array<IndexRange,5>&, const double th = 1.0e-8);

    // add RDMs (using fock-multiplied 4RDM)
    void append(const int jst, const int ist, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                              std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    // diagonalize and set to shalf and denom
    void compute();

    const VecView denom_x() const { return denom_x_; }
    const VecView denom_h() const { return denom_h_; }
    const VecView denom_xx() const { return denom_xx_; }
    const VecView denom_hh() const { return denom_hh_; }
    const VecView denom_xh() const { return denom_xh_; }
    const VecView denom_xhh() const { return denom_xhh_; }
    const VecView denom_xxh() const { return denom_xxh_; }

    std::shared_ptr<const TATensor<DataType,2>> tashalf_x() const { return tashalf_x_; }
    std::shared_ptr<const TATensor<DataType,2>> tashalf_h() const { return tashalf_h_; }
    std::shared_ptr<const TATensor<DataType,3>> tashalf_xx() const { return tashalf_xx_; }
    std::shared_ptr<const TATensor<DataType,3>> tashalf_hh() const { return tashalf_hh_; }
    std::shared_ptr<const TATensor<DataType,4>> tashalf_xhh() const { return tashalf_xhh_; }
    std::shared_ptr<const TATensor<DataType,4>> tashalf_xxh() const { return tashalf_xxh_; }
    template<bool I>
    std::shared_ptr<const TATensor<DataType,3>> tashalf_xh() const { auto o = I ? tashalf_xh_ : tashalf_xh2_; assert(o); return o; }
};

extern template class Denom<double>;
extern template class Denom<std::complex<double>>;

}
}

#endif

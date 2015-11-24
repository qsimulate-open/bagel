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
    const double thresh_;

    std::shared_ptr<MatType> shalf_x_;
    std::shared_ptr<MatType> shalf_h_;
    std::shared_ptr<MatType> shalf_xx_;
    std::shared_ptr<MatType> shalf_hh_;
    std::shared_ptr<MatType> shalf_xh_;
    std::shared_ptr<MatType> shalf_xhh_;
    std::shared_ptr<MatType> shalf_xxh_;

    std::shared_ptr<MatType> work_x_;
    std::shared_ptr<MatType> work_h_;
    std::shared_ptr<MatType> work_xx_;
    std::shared_ptr<MatType> work_hh_;
    std::shared_ptr<MatType> work_xh_;
    std::shared_ptr<MatType> work_xhh_;
    std::shared_ptr<MatType> work_xxh_;

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
    Denom(std::shared_ptr<const MatType> fock, const int nstates, const double th = 1.0e-8);

    // add RDMs (using fock-multiplied 4RDM)
    void append(const int jst, const int ist, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                              std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    // add RDMs (using original 4RDM)
    void append(const int jst, const int ist, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                              std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<4,DataType>>);
    // add RDMs (using Kramers-reduced 4RDM)
    void append(const int jst, const int ist, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                              std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const Kramers<8,RDM<4,DataType>>>);
    // diagonalize and set to shalf and denom
    void compute();

    const VectorB& denom_x() const { return denom_x_; }
    const VectorB& denom_h() const { return denom_h_; }
    const VectorB& denom_xx() const { return denom_xx_; }
    const VectorB& denom_hh() const { return denom_hh_; }
    const VectorB& denom_xh() const { return denom_xh_; }
    const VectorB& denom_xhh() const { return denom_xhh_; }
    const VectorB& denom_xxh() const { return denom_xxh_; }

    std::shared_ptr<const TATensor<DataType,4>> tashalf_xxh(const std::vector<IndexRange>& ranges) const {
      auto out = std::make_shared<TATensor<DataType,4>>(ranges);
      // TODO not an optimal code
      auto tmp = std::make_shared<btas::Tensor4<DataType>>(btas::CRange<3>(ranges[0].size(), ranges[1].size(), ranges[2].size(), ranges[3].size()));
      std::copy_n(shalf_xxh_->data(), shalf_xxh_->size(), tmp->data());
      const int nclo = ranges[1].front().offset();
      fill_block<4,DataType>(out, tmp, std::vector<int>{nclo, nclo, nclo, 0});
      return out;
    }
    std::shared_ptr<const TATensor<DataType,4>> tashalf_xhh(const std::vector<IndexRange>& ranges) const {
      auto out = std::make_shared<TATensor<DataType,4>>(ranges);
      // TODO not an optimal code
      auto tmp = std::make_shared<btas::Tensor4<DataType>>(btas::CRange<3>(ranges[0].size(), ranges[1].size(), ranges[2].size(), ranges[3].size()));
      std::copy_n(shalf_xhh_->data(), shalf_xhh_->size(), tmp->data());
      const int nclo = ranges[1].front().offset();
      fill_block<4,DataType>(out, tmp, std::vector<int>{nclo, nclo, nclo, 0});
      return out;
    }
    std::shared_ptr<const TATensor<DataType,3>> tashalf_xx(const std::vector<IndexRange>& ranges) const {
      auto out = std::make_shared<TATensor<DataType,3>>(ranges);
      // TODO not an optimal code
      auto tmp = std::make_shared<btas::Tensor3<DataType>>(btas::CRange<3>(ranges[0].size(), ranges[1].size(), ranges[2].size()));
      std::copy_n(shalf_xx_->data(), shalf_xx_->size(), tmp->data());
      const int nclo = ranges[1].front().offset();
      fill_block<3,DataType>(out, tmp, std::vector<int>{nclo, nclo, 0});
      return out;
    }
    std::shared_ptr<const TATensor<DataType,3>> tashalf_hh(const std::vector<IndexRange>& ranges) const {
      auto out = std::make_shared<TATensor<DataType,3>>(ranges);
      // TODO not an optimal code
      auto tmp = std::make_shared<btas::Tensor3<DataType>>(btas::CRange<3>(ranges[0].size(), ranges[1].size(), ranges[2].size()));
      std::copy_n(shalf_hh_->data(), shalf_hh_->size(), tmp->data());
      const int nclo = ranges[1].front().offset();
      fill_block<3,DataType>(out, tmp, std::vector<int>{nclo, nclo, 0});
      return out;
    }
    std::shared_ptr<const TATensor<DataType,2>> tashalf_x(const std::vector<IndexRange>& ranges) const {
      auto out = std::make_shared<TATensor<DataType,2>>(ranges);
      const int nclo = ranges[1].front().offset();
      fill_block<2,DataType>(out, shalf_x_, std::vector<int>{nclo, 0});
      return out;
    }
    std::shared_ptr<const TATensor<DataType,2>> tashalf_h(const std::vector<IndexRange>& ranges) const {
      auto out = std::make_shared<TATensor<DataType,2>>(ranges);
      const int nclo = ranges[1].front().offset();
      fill_block<2,DataType>(out, shalf_h_, std::vector<int>{nclo, 0});
      return out;
    }

    // deprecated
    std::shared_ptr<const MatType> shalf_x() const { return shalf_x_; }
    std::shared_ptr<const MatType> shalf_h() const { return shalf_h_; }
    std::shared_ptr<const MatType> shalf_xx() const { return shalf_xx_; }
    std::shared_ptr<const MatType> shalf_hh() const { return shalf_hh_; }
    std::shared_ptr<const MatType> shalf_xh() const { return shalf_xh_; }
    std::shared_ptr<const MatType> shalf_xhh() const { return shalf_xhh_; }
    std::shared_ptr<const MatType> shalf_xxh() const { return shalf_xxh_; }
    // deprecated
    const double& denom_x(const size_t i) const { return denom_x_(i); }
    const double& denom_h(const size_t i) const { return denom_h_(i); }
    const double& denom_xx(const size_t i) const { return denom_xx_(i); }
    const double& denom_hh(const size_t i) const { return denom_hh_(i); }
    const double& denom_xh(const size_t i) const { return denom_xh_(i); }
    const double& denom_xhh(const size_t i) const { return denom_xhh_(i); }
    const double& denom_xxh(const size_t i) const { return denom_xxh_(i); }
};

template<>
void Denom<double>::append(const int, const int, std::shared_ptr<const RDM<1>>, std::shared_ptr<const RDM<2>>,
                                                 std::shared_ptr<const RDM<3>>, std::shared_ptr<const Kramers<8,RDM<4>>>);
template<>
void Denom<std::complex<double>>::append(const int, const int, std::shared_ptr<const ZRDM<1>>, std::shared_ptr<const ZRDM<2>>,
                                                               std::shared_ptr<const ZRDM<3>>, std::shared_ptr<const Kramers<8,ZRDM<4>>>);

extern template class Denom<double>;
extern template class Denom<std::complex<double>>;

}
}

#endif

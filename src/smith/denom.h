//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: denom.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_SMITH_DENOM_H
#define __SRC_SMITH_DENOM_H

#include <src/wfn/rdm.h>
#include <src/util/math/zmatrix.h>

namespace bagel {
namespace SMITH {

template<typename DataType>
class Denom {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
    using ViewType = typename std::conditional<std::is_same<DataType,double>::value,MatView,ZMatView>::type;

  protected:
    std::shared_ptr<const MatType> fock_;
    const double thresh_;
    const size_t nstates_;
    std::map<std::string,size_t> dim_;

    // init functions
    void init_(const std::string& tag, const int, const int,
               std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
               std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    virtual void set_sblock(const std::string tag, const int, const int, const btas::TensorView2<DataType> data) = 0;
    virtual void set_wblock(const std::string tag, const int, const int, const btas::TensorView2<DataType> data) = 0;

  public:
    Denom(std::shared_ptr<const MatType> fock, const int nstates, const double thresh_overlap);

    // add RDMs (using fock-multiplied 4RDM)
    void append(const int jst, const int ist, std::shared_ptr<const RDM<1,DataType>>, std::shared_ptr<const RDM<2,DataType>>,
                                              std::shared_ptr<const RDM<3,DataType>>, std::shared_ptr<const RDM<3,DataType>>);
    // diagonalize and set to shalf and denom
    virtual void compute() = 0;

    virtual const ViewType shalf(const std::string, const int) const = 0;
    virtual double denom(const std::string& tag, const int ist, const size_t i) const = 0;
};


template<typename DataType>
class Denom_SSSR : public Denom<DataType> {
  protected:
    using Denom<DataType>::thresh_;
    using Denom<DataType>::nstates_;
    using Denom<DataType>::dim_;
    using Denom<DataType>::fock_;
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
    using ViewType = typename std::conditional<std::is_same<DataType,double>::value,MatView,ZMatView>::type;

  protected:
    std::map<std::string,std::vector<std::shared_ptr<MatType>>> shalf_;
    std::map<std::string,std::vector<std::shared_ptr<MatType>>> work_;
    std::map<std::string,std::vector<VectorB>> denom_;

    void set_sblock(const std::string tag, const int jst, const int ist, const btas::TensorView2<DataType> o) override {
      assert(jst == ist && o.size() == shalf_.at(tag)[ist]->size());
      std::copy_n(&*o.begin(), o.size(), shalf_.at(tag)[ist]->data());
    }
    void set_wblock(const std::string tag, const int jst, const int ist, const btas::TensorView2<DataType> o) override {
      assert(jst == ist && o.size() == work_.at(tag)[ist]->size());
      std::copy_n(&*o.begin(), o.size(), work_.at(tag)[ist]->data());
    }

  public:
    Denom_SSSR(std::shared_ptr<const MatType> fock, const int nstates, const double thresh_overlap);

    const ViewType shalf(const std::string tag, const int ist) const { return *(shalf_.at(tag)[ist]); }
    void compute() override;
    double denom(const std::string& tag, const int ist, const size_t i) const override { return denom_.at(tag)[ist](i); }
};


template<typename DataType>
class Denom_MSMR : public Denom<DataType> {
  protected:
    using Denom<DataType>::thresh_;
    using Denom<DataType>::dim_;
    using Denom<DataType>::fock_;
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;
    using ViewType = typename std::conditional<std::is_same<DataType,double>::value,MatView,ZMatView>::type;

  protected:
    std::map<std::string,std::shared_ptr<MatType>> shalf_;
    std::map<std::string,std::shared_ptr<MatType>> work_;
    std::map<std::string, VectorB> denom_;

    void set_sblock(const std::string tag, const int jst, const int ist, const btas::TensorView2<DataType> o) override {
      const size_t dim = dim_.at(tag);
      shalf_.at(tag)->copy_block(jst*dim, ist*dim, dim, dim, o);
    }
    void set_wblock(const std::string tag, const int jst, const int ist, const btas::TensorView2<DataType> o) override {
      const size_t dim = dim_.at(tag);
      work_.at(tag)->copy_block(jst*dim, ist*dim, dim, dim, o);
    }

  public:
    Denom_MSMR(std::shared_ptr<const MatType> fock, const int nstates, const double thresh_overlap);

    const ViewType shalf(const std::string tag, const int ist) const {
      const size_t dim = dim_.at(tag);
      return shalf_.at(tag)->slice(ist*dim, (ist+1)*dim);
    }
    void compute() override;
    double denom(const std::string& tag, const int ist, const size_t i) const override { return denom_.at(tag)(i); }
};

extern template class Denom<double>;
extern template class Denom<std::complex<double>>;

}
}

#endif

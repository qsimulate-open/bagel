//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moint.h
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

// A base class for electorn Correlation methods
// Certain kinds of MO integrals are formed.
//   - aaii (assumes DF - TODO half transformed DF vector might be available..)
//

#ifndef __SRC_SMITH_MOINT_H
#define __SRC_SMITH_MOINT_H

#include <stddef.h>
#include <memory>
#include <stdexcept>
#include <src/smith/smith_info.h>
#include <src/smith/tensor.h>
#include <src/scf/hf/fock.h>
#include <src/util/math/algo.h>

namespace bagel {
namespace SMITH {

// the template parameter T specifies the storage type

template<typename DataType>
class K2ext {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;

    std::shared_ptr<const SMITH_Info<DataType>> info_;
    std::shared_ptr<const MatType> coeff_;
    std::vector<IndexRange> blocks_;
    std::shared_ptr<Tensor_<DataType>> data_;

    // some handwritten drivers
    void init() { assert(false); }

  public:
    K2ext(std::shared_ptr<const SMITH_Info<DataType>> r, std::shared_ptr<const MatType> c, const std::vector<IndexRange>& b);

    std::shared_ptr<Tensor_<DataType>> tensor() { return data_; }
};
template<> void K2ext<double>::init();
template<> void K2ext<std::complex<double>>::init();
extern template class K2ext<double>;
extern template class K2ext<std::complex<double>>;


template<typename DataType>
class MOFock {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value,Matrix,ZMatrix>::type;

    std::shared_ptr<const SMITH_Info<DataType>> info_;
    std::shared_ptr<const MatType> coeff_;
    std::vector<IndexRange> blocks_;
    std::shared_ptr<Tensor_<DataType>> data_;
    std::shared_ptr<Tensor_<DataType>> h1_;

    double core_energy_;

    void init() { assert(false); }

  public:
    MOFock(std::shared_ptr<const SMITH_Info<DataType>> r, const std::vector<IndexRange>& b);

    // fock operator
    std::shared_ptr<Tensor_<DataType>> tensor() { return data_; }
    // core Fock operator minus diagonal part of the two-body integrals
    std::shared_ptr<Tensor_<DataType>> h1() { return h1_; }

    std::shared_ptr<const MatType> coeff() const { return coeff_; }
    double core_energy() const { return core_energy_; }
};
template<> void MOFock<double>::init();
template<> void MOFock<std::complex<double>>::init();
extern template class MOFock<double>;
extern template class MOFock<std::complex<double>>;


}
}

#endif


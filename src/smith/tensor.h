//
// BAGEL - Parallel electron correlation program.
// Filename: tensor.h
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


#ifndef __SRC_SMITH_TENSOR_H
#define __SRC_SMITH_TENSOR_H

#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <unordered_set>
#include <src/ci/fci/civec.h>
#include <src/util/math/matrix.h>
#include <src/util/math/matop.h>
#include <src/util/prim_op.h>
#include <src/smith/storage.h>
#include <src/smith/indexrange.h>
#include <src/smith/loopgenerator.h>

namespace bagel {
namespace SMITH {

template <typename DataType>
class Tensor_ {
  protected:
    using MatType = typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type;
  protected:
    std::vector<IndexRange> range_;
    std::shared_ptr<Storage<DataType>> data_;
    int rank_;
    std::unordered_set<size_t> sparse_;

    mutable bool initialized_;

    bool allocated_;

  public:
    Tensor_(std::vector<IndexRange> in, const bool kramers = false, const std::unordered_set<size_t> sparse = {}, const bool alloc = false);

    Tensor_<DataType>& operator=(const Tensor_<DataType>& o) {
      *data_ = *(o.data_);
      allocated_ = true;
      return *this;
    }

    std::shared_ptr<Tensor_<DataType>> clone() const {
      return std::make_shared<Tensor_<DataType>>(range_, false, sparse_);
    }

    std::shared_ptr<Tensor_<DataType>> copy() const {
      std::shared_ptr<Tensor_<DataType>> out = clone();
      *out = *this;
      return out;
    }

    virtual void init() const { initialized_ = true; }

    void ax_plus_y(const DataType& a, const Tensor_<DataType>& o) { data_->ax_plus_y(a, o.data_); }
    void ax_plus_y(const DataType& a, std::shared_ptr<const Tensor_<DataType>> o) { ax_plus_y(a, *o); }

    void scale(const DataType& a) { data_->scale(a); }

    bool allocated() const { return allocated_; }
    void allocate();

    template<typename ...args>
    bool is_local(args&& ...p) const { return data_->is_local(std::forward<args>(p)...); }
    template<typename ...args>
    bool exists(args&& ...p) const { return sparse_.empty() || sparse_.count(generate_hash_key(std::forward<args>(p)...)); }

    DataType dot_product(const Tensor_<DataType>& o) const { return data_->dot_product(*o.data_); }
    DataType dot_product(std::shared_ptr<const Tensor_<DataType>> o) const { return dot_product(*o); }

    int rank() const { return rank_; }
    size_t size_alloc() const;

    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double rms() const { return std::sqrt(detail::real(dot_product(*this))/size_alloc()); }

    std::vector<IndexRange> indexrange() const { return range_; }

    template<typename ...args>
    std::unique_ptr<DataType[]> get_block(args&& ...p) const {
      return data_->get_block(std::forward<args>(p)...);
    }

    template<typename ...args>
    void put_block(std::unique_ptr<DataType[]>& o, args&& ...p) {
      data_->put_block(o, std::forward<args>(p)...);
    }

    template<typename ...args>
    void add_block(std::unique_ptr<DataType[]>& o, args&& ...p) {
      data_->add_block(o, std::forward<args>(p)...);
    }

    template<typename ...args>
    size_t get_size(args&& ...p) const {
      return data_->blocksize(std::forward<args>(p)...);
    }

    void zero() {
      data_->zero();
    }

    std::vector<DataType> diag() const;

    std::shared_ptr<MatType> matrix() const;
    std::shared_ptr<MatType> matrix2() const;

    std::shared_ptr<Civector<DataType>> civec(std::shared_ptr<const Determinants> det) const;

    // for Kramers tensors (does not do anything for standard tensors)
    void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& p) { data_->set_perm(p); }
    void set_stored_sectors(const std::list<std::vector<bool>>& s) { data_->set_stored_sectors(s); }
};

extern template class Tensor_<double>;
extern template class Tensor_<std::complex<double>>;

namespace CASPT2 { using Tensor = Tensor_<double>; }
namespace MRCI   { using Tensor = Tensor_<double>; }
namespace RelCASPT2 { using Tensor = Tensor_<std::complex<double>>; }
namespace RelMRCI   { using Tensor = Tensor_<std::complex<double>>; }

}
}

#endif

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
#include <list>
#include <map>
#include <memory>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <type_traits>
#include <src/ci/fci/civec.h>
#include <src/util/math/matrix.h>
#include <src/util/math/matop.h>
#include <src/util/prim_op.h>
#include <src/smith/storage.h>
#include <src/smith/indexrange.h>
#include <src/smith/loopgenerator.h>
#include <src/smith/tatensor.h>

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

    virtual void init() const { initialized_ = true; }
    mutable bool initialized_;

  public:
    Tensor_(std::vector<IndexRange> in, const bool kramers = false);

    template<typename D, int N, class = typename std::enable_if<std::is_same<D,DataType>::value>::type>
    Tensor_(const TATensor<D,N>& o) : Tensor_(o.indexrange()) { // delegate constuctor
      // loop over tiles
      for (auto it = o.begin(); it != o.end(); ++it) {
        // first get range
        const TiledArray::Range range = o.trange().make_tile_range(it.ordinal());
        auto lo = range.lobound();
        size_t cnt = 0;
        std::vector<Index> indices;
        for (auto i = lo.rbegin(); i != lo.rend(); ++i, ++cnt) {
          Index ind;
          for (auto& index : range_[cnt])
            if (index.offset() == *i) {
              ind = index;
              break;
            }
          assert(ind.size());
          indices.push_back(ind);
        }
        const typename TiledArray::Array<DataType,N>::value_type& tile = *it;
        std::unique_ptr<DataType[]> data(new DataType[tile.size()]);
        std::copy_n(&(tile[0]), tile.size(), data.get());
        this->put_block(data, indices);
      }
    }

    Tensor_<DataType>& operator=(const Tensor_<DataType>& o) {
      *data_ = *(o.data_);
      return *this;
    }

    // Debug function to return TiledArray Tensor from this
    template<int N>
    std::shared_ptr<TATensor<DataType,N>> tiledarray() {
      assert(range_.size() == N);
      std::vector<std::map<size_t,size_t>> keymap;
      for (auto it = range_.rbegin(); it != range_.rend(); ++it) {
        std::map<size_t,size_t> key;
        for (auto& j : *it)
          key.emplace(j.offset(), j.key());
        keymap.push_back(key);
      }
      auto out = std::make_shared<TATensor<DataType,N>>(range_);

      for (auto it = out->begin(); it != out->end(); ++it) {
        const TiledArray::Range range = out->trange().make_tile_range(it.ordinal());
        typename TiledArray::Array<DataType,N>::value_type tile(range);
        auto lo = range.lobound();
        assert(lo.size() == N);
        std::vector<size_t> seed(lo.size());
        {
          auto iter = seed.rbegin();
          auto key = keymap.begin();
          for (auto jt = lo.begin(); jt != lo.end(); ++jt, ++key, ++iter)
            *iter = key->at(*jt);
        }
        // pull out the tile
        std::unique_ptr<DataType[]> data = move_block(generate_hash_key(seed));
        // copy
        std::copy_n(data.get(), tile.size(), &(tile[0]));
        *it = tile;
        // restore the original tensnor
        put_block(data, generate_hash_key(seed));
      }

      return out;
    }

    std::shared_ptr<Tensor_<DataType>> clone() const {
      return std::make_shared<Tensor_<DataType>>(range_);
    }

    std::shared_ptr<Tensor_<DataType>> copy() const {
      std::shared_ptr<Tensor_<DataType>> out = clone();
      *out = *this;
      return out;
    }

    void ax_plus_y(const DataType& a, const Tensor_<DataType>& o) { data_->ax_plus_y(a, o.data_); }
    void ax_plus_y(const DataType& a, std::shared_ptr<const Tensor_<DataType>> o) { ax_plus_y(a, *o); }

    void scale(const DataType& a) { data_->scale(a); }

    DataType dot_product(const Tensor_<DataType>& o) const { return data_->dot_product(*o.data_); }
    DataType dot_product(std::shared_ptr<const Tensor_<DataType>> o) const { return dot_product(*o); }

    int rank() const { return rank_; }
    size_t size_alloc() const;

    double norm() const { return std::sqrt(detail::real(dot_product(*this))); }
    double rms() const { return std::sqrt(detail::real(dot_product(*this))/size_alloc()); }

    std::vector<IndexRange> indexrange() const { return range_; }

    template<typename ...args>
    std::unique_ptr<DataType[]> get_block(args&& ...p) const {
      if (!initialized_) init();
      return data_->get_block(std::forward<args>(p)...);
    }

    template<typename ...args>
    std::unique_ptr<DataType[]> move_block(args&& ...p) {
      return data_->move_block(std::forward<args>(p)...);
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

    template<typename ...args>
    size_t get_size_alloc(args&& ...p) const {
      return data_->blocksize_alloc(std::forward<args>(p)...);
    }

    void zero() {
      data_->zero();
    }

    void conjugate_inplace() {
      data_->conjugate_inplace();
    }

    std::vector<DataType> diag() const;

    std::shared_ptr<MatType> matrix() const;
    std::shared_ptr<MatType> matrix2() const;

    std::shared_ptr<Civector<DataType>> civec(std::shared_ptr<const Determinants> det) const;

    // for Kramers tensors (does not do anything for standard tensors)
    void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& p) { data_->set_perm(p); }

    void print1(std::string label, const double thresh = 5.0e-2) const;
    void print2(std::string label, const double thresh = 5.0e-2) const;
    void print3(std::string label, const double thresh = 5.0e-2) const;
    void print4(std::string label, const double thresh = 5.0e-2) const;
    void print5(std::string label, const double thresh = 5.0e-2) const;
    void print6(std::string label, const double thresh = 5.0e-2) const;
    void print8(std::string label, const double thresh = 5.0e-2) const;
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

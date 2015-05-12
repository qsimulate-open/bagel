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

namespace bagel {
namespace SMITH {

// this assumes < 256 blocks; TODO runtime determination?
const static int shift = 8;

static size_t generate_hash_key(const std::vector<size_t>& o) {
  size_t out = 0;
  for (auto i = o.rbegin(); i != o.rend(); ++i) { out <<= shift; out += *i; }
  return out;
}

static size_t generate_hash_key() { return 0; }

template<class T, typename... args>
size_t generate_hash_key(const T& head, const args&... tail) {
  return (generate_hash_key(tail...) << shift) + head.key();
}


class Tensor {
  protected:

    std::vector<IndexRange> range_;
    std::shared_ptr<Storage> data_;
    int rank_;

    virtual void init() const { initialized_ = true; }
    mutable bool initialized_;

  public:
    Tensor(std::vector<IndexRange> in);

    Tensor& operator=(const Tensor& o) {
      *data_ = *(o.data_);
      return *this;
    }

    std::shared_ptr<Tensor> clone() const {
      return std::make_shared<Tensor>(range_);
    }

    std::shared_ptr<Tensor> copy() const {
      std::shared_ptr<Tensor> out = clone();
      *out = *this;
      return out;
    }

    void ax_plus_y(const double a, const Tensor& o) { data_->ax_plus_y(a, o.data_); }
    void ax_plus_y(const double a, const std::shared_ptr<Tensor> o) { data_->ax_plus_y(a, o->data_); }

    void scale(const double a) { data_->scale(a); }

    double dot_product(const Tensor& o) const { return data_->dot_product(*o.data_); }
    double dot_product(const std::shared_ptr<Tensor>& o) const { return data_->dot_product(*o->data_); }

    size_t size_alloc() const;

    double norm() const { return std::sqrt(dot_product(*this)); }
    double rms() const { return std::sqrt(dot_product(*this)/size_alloc()); }

    std::vector<IndexRange> indexrange() const { return range_; }

    template<typename ...args>
    std::unique_ptr<double[]> get_block(const args& ...p) const {
      if (!initialized_) init();
      return data_->get_block(generate_hash_key(p...));
    }

    template<typename ...args>
    std::unique_ptr<double[]> move_block(const args& ...p) {
      return data_->move_block(generate_hash_key(p...));
    }

    template<typename ...args>
    void put_block(std::unique_ptr<double[]>& o, const args& ...p) {
      data_->put_block(generate_hash_key(p...), o);
    }

    template<typename ...args>
    void add_block(std::unique_ptr<double[]>& o, const args& ...p) {
      data_->add_block(generate_hash_key(p...), o);
    }

    template<typename ...args>
    size_t get_size(const args& ...p) const {
      return data_->blocksize(generate_hash_key(p...));
    }

    template<typename ...args>
    size_t get_size_alloc(const args& ...p) const {
      return data_->blocksize_alloc(generate_hash_key(p...));
    }

    void zero() {
      data_->zero();
    }

    std::vector<double> diag() const;

    std::shared_ptr<Matrix> matrix() const;
    std::shared_ptr<Matrix> matrix2() const;


    std::shared_ptr<Civec> civec(std::shared_ptr<const Determinants> det) const;


    void print1(std::string label, const double thresh = 5.0e-2) const;
    void print2(std::string label, const double thresh = 5.0e-2) const;
    void print3(std::string label, const double thresh = 5.0e-2) const;
    void print4(std::string label, const double thresh = 5.0e-2) const;
    void print5(std::string label, const double thresh = 5.0e-2) const;
    void print6(std::string label, const double thresh = 5.0e-2) const;
    void print8(std::string label, const double thresh = 5.0e-2) const;
};

}
}

#endif

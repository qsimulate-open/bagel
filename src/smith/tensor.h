//
// Newint - Parallel electron correlation program.
// Filename: tensor.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_SMITH_TENSOR_H
#define __SRC_SMITH_TENSOR_H

#include <list>
#include <map>
#include <memory>
#include <iostream>
#include <cassert>
#include <src/smith/storage.h>
#include <src/smith/indexrange.h>
#include <src/smith/loopgenerator.h>

namespace SMITH {

// this function should be fast
inline static size_t generate_hash_key(const std::vector<size_t>& o) {
  // this assumes < 256 blocks; TODO runtime determination?
  const int shift = 8;
  size_t out = 0;
  for (auto i = o.begin(); i != o.end(); ++i) {
    out <<= shift;
    out ^= (*i);
    assert(*i < (1 << shift));
  }
  return out;
};

template <typename T>
class Tensor {
  protected:
    std::vector<IndexRange> range_;
    std::shared_ptr<T> data_; 
    int rank_;
    bool initialized_;

  public:
    Tensor(std::vector<IndexRange> in, bool init = true) : range_(in), rank_(in.size()), initialized_(init) {
      // make blocl list
      if (!in.empty()) {
        LoopGenerator lg(in);
        std::vector<std::vector<Index> > index = lg.block_loop();

        // first compute hashtags and length 
        std::map<size_t, size_t> hashmap;
        size_t off = 0;
        for (auto i = index.begin(); i != index.end(); ++i) {
          size_t size = 1lu;
          std::vector<size_t> h;
          for (auto j = i->begin(); j != i->end(); ++j) {
            size *= j->size();
            h.push_back(j->key());
          }
          hashmap.insert(std::make_pair(generate_hash_key(h), size)); 
          off += size;
        }

        std::shared_ptr<T> tmp(new T(hashmap, init));
        data_ = tmp;
      } else {
        rank_ = 1;
        std::map<size_t, size_t> hashmap;
        hashmap.insert(std::make_pair(0lu, 1lu));
        std::shared_ptr<T> tmp(new T(hashmap, init));
        data_ = tmp;
      }
    };

    void initialize() { data_->initialize(); };

    ~Tensor() {};

    Tensor<T>& operator=(const Tensor<T>& o) {
      *data_ = *(o.data_);
      return *this;
    };

    std::shared_ptr<Tensor<T> > clone() const {
      std::shared_ptr<Tensor<T> > out(new Tensor<T>(range_));
      return out;
    };
    std::shared_ptr<Tensor<T> > copy() const {
      std::shared_ptr<Tensor<T> > out = clone();
      *out = *this;
      return out;
    };

    void daxpy(const double a, const Tensor<T>& o) { data_->daxpy(a, o.data_); };
    void daxpy(const double a, const std::shared_ptr<Tensor<T> > o) { data_->daxpy(a, o->data_); };

    void scale(const double a) { data_->scale(a); };

    double ddot(const Tensor<T>& o) { return data_->ddot(*o.data_); };
    double ddot(const std::shared_ptr<Tensor<T> >& o) { return data_->ddot(*o->data_); };

    size_t size() const { return data_->length(); };
    size_t length() const { return data_->length(); };

    double norm() { return std::sqrt(ddot(*this)); };
    double rms() { return std::sqrt(ddot(*this)/size()); };

    std::vector<IndexRange> indexrange() const { return range_; };

    std::unique_ptr<double[]> get_block(const std::vector<size_t>& p) const {
      assert(p.size() == rank_); 
      if (!data_) throw std::logic_error("Tensor not initialized");
      return std::move(data_->get_block(generate_hash_key(p)));
    };

    std::unique_ptr<double[]> move_block(const std::vector<size_t>& p) {
      assert(p.size() == rank_);
      return std::move(data_->move_block(generate_hash_key(p)));
    };

    void put_block(const std::vector<size_t>& p, std::unique_ptr<double[]>& o) {
      data_->put_block(generate_hash_key(p), o);
    };

    void add_block(const std::vector<size_t>& p, const std::unique_ptr<double[]>& o) {
      if (!data_) throw std::logic_error("Tensor not initialized");
      data_->add_block(generate_hash_key(p), o);
    };

    size_t get_size(const std::vector<size_t>& p) {
      assert(p.size() == rank_); 
      return data_->blocksize(generate_hash_key(p));
    };

    void zero() {
      data_->zero();
    };

    std::unique_ptr<double[]> diag() {
      if (rank_ != 2 || range_[0] != range_[1])
        throw std::logic_error("Tensor::diag can be called only with a square tensor of rank 2");
      const size_t size = range_[0].size();
      std::unique_ptr<double[]> buf(new double[size]);
      for (auto i = range_[0].begin(); i != range_[0].end(); ++i) {
        std::vector<size_t> h = vec(i->key(), i->key());
        std::unique_ptr<double[]> data0 = move_block(h); 
        for (int j = 0; j != i->size(); ++j) {
          buf[j+i->offset()] = data0[j+j*i->size()]; 
        }
        put_block(h, data0);
      }
      return std::move(buf);
    };

    std::shared_ptr<Tensor<T> > add_dagger() {
      std::shared_ptr<Tensor<T> > out = clone();
      std::vector<IndexRange> o = indexrange();
      assert(o.size() == 4);
      for (auto i3 = o[3].range().begin(); i3 != o[3].range().end(); ++i3) {
        for (auto i2 = o[2].range().begin(); i2 != o[2].range().end(); ++i2) {
          for (auto i1 = o[1].range().begin(); i1 != o[1].range().end(); ++i1) {
            for (auto i0 = o[0].range().begin(); i0 != o[0].range().end(); ++i0) {
              std::vector<size_t> h(4); h[0] = i0->key(); h[1] = i1->key(); h[2] = i2->key(); h[3] = i3->key();
              std::vector<size_t> g(4); g[0] = i2->key(); g[1] = i3->key(); g[2] = i0->key(); g[3] = i1->key();
              std::unique_ptr<double[]> data0 = get_block(h);
              const std::unique_ptr<double[]> data1 = get_block(g);
              sort_indices<2,3,0,1,1,1,1,1>(data1, data0, i2->size(), i3->size(), i0->size(), i1->size()); 
              out->put_block(h,data0);
            }
          }
        }
      }
      return out;
    };
};

}

#endif


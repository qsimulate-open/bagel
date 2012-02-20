//
// Newint - Parallel electron correlation program.
// Filename: tensor.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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
    const int rank_;

  public:
    Tensor(std::vector<IndexRange> in) : range_(in), rank_(in.size()) {
      // make blocl list
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

      std::shared_ptr<T> tmp(new T(hashmap));
      data_ = tmp;
    };

    ~Tensor() {};

    Tensor<T>& operator=(const Tensor<T>& o) {
      *data_ = *(o.data_);
      return *this;
    };

    std::shared_ptr<Tensor<T> > clone() const {
      std::shared_ptr<Tensor<T> > out(new Tensor<T>(range_));
      return out;
    };

    std::vector<IndexRange> indexrange() const { return range_; };

    std::unique_ptr<double[]> get_block(const std::vector<size_t>& p) const {
      assert(p.size() == rank_); 
      return std::move(data_->get_block(generate_hash_key(p)));
    };

    void put_block(const std::vector<size_t>& p, std::unique_ptr<double[]>& o) {
      data_->put_block(generate_hash_key(p), o);
    };

    void add_block(const std::vector<size_t>& p, const std::unique_ptr<double[]>& o) {
      data_->add_block(generate_hash_key(p), o);
    };

    size_t get_size(const std::vector<size_t>& p) {
      assert(p.size() == rank_); 
      return data_->blocksize(generate_hash_key(p));
    };

    void zero() { data_->zero(); };
};

}

#endif


//
// Author : Toru Shiozaki
// Date   : Feb 2012
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
    std::shared_ptr<T> data_; 
    const int rank_;

  public:
    Tensor(std::vector<IndexRange> in) : rank_(in.size()) {
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


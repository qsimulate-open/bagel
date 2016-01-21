//
// BAGEL - Parallel electron correlation program.
// Filename: storage.h
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

//
// this is a base class of storage object
// ... It can be distributed memory, distributed disk, etc.
// ... All of them should be derived from here with the same interface.

#ifndef __SRC_SMITH_STORAGE_H
#define __SRC_SMITH_STORAGE_H

#include <stddef.h>
#include <map>
#include <unordered_map>
#include <memory>
#include <tuple>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <complex>
#include <src/smith/indexrange.h>

namespace bagel {
namespace SMITH {


// this assumes < 256 blocks; TODO runtime determination?
const static int shift = 8;

static size_t generate_hash_key(const std::vector<size_t>& o) {
  size_t out = 0;
  for (auto i = o.rbegin(); i != o.rend(); ++i) { out <<= shift; out += *i; }
  return out;
}

static size_t generate_hash_key(const std::vector<Index>& o) {
  size_t out = 0;
  for (auto i = o.rbegin(); i != o.rend(); ++i) { out <<= shift; out += i->key(); }
  return out;
}

static size_t generate_hash_key() { return 0; }

template<class T, typename... args>
size_t generate_hash_key(const T& head, const args&... tail) {
  return (generate_hash_key(tail...) << shift) + head.key();
}


template<typename DataType>
class StorageIncore {
  protected:
    // Global Array handler:
    int ga_;
    int64_t totalsize_;

    // this relates hash keys to lo and high of the block
    std::unordered_map<size_t, std::pair<int64_t, int64_t>> hashtable_;
    // distribution information
    std::vector<int64_t> blocks_;

    bool initialized_;

    std::unique_ptr<DataType[]> get_block_(const size_t& key) const;
    std::unique_ptr<DataType[]> move_block_(const size_t& key);
    void put_block_(std::unique_ptr<DataType[]>& dat, const size_t& key);
    void add_block_(const std::unique_ptr<DataType[]>& dat, const size_t& key);

  public:
    StorageIncore(const std::map<size_t, size_t>& size, bool init);
    ~StorageIncore();

    virtual std::unique_ptr<DataType[]> get_block() const;
    virtual std::unique_ptr<DataType[]> get_block(const Index& i0) const;
    virtual std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1) const;
    virtual std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2) const;
    virtual std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3) const;
    virtual std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                  const Index& i4) const;
    virtual std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                  const Index& i4, const Index& i5) const;
    virtual std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                  const Index& i4, const Index& i5, const Index& i6) const;
    virtual std::unique_ptr<DataType[]> get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                  const Index& i4, const Index& i5, const Index& i6, const Index& i7) const;
    virtual std::unique_ptr<DataType[]> get_block(std::vector<Index> i) const;

    virtual std::unique_ptr<DataType[]> move_block();
    virtual std::unique_ptr<DataType[]> move_block(const Index& i0);
    virtual std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1);
    virtual std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2);
    virtual std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3);
    virtual std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                   const Index& i4);
    virtual std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                   const Index& i4, const Index& i5);
    virtual std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                   const Index& i4, const Index& i5, const Index& i6);
    virtual std::unique_ptr<DataType[]> move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                   const Index& i4, const Index& i5, const Index& i6, const Index& i7);

    virtual void put_block(std::unique_ptr<DataType[]>& dat);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                             const Index& i4);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                             const Index& i4, const Index& i5);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                             const Index& i4, const Index& i5, const Index& i6);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                             const Index& i4, const Index& i5, const Index& i6, const Index& i7);
    virtual void put_block(std::unique_ptr<DataType[]>& dat, std::vector<Index> i);

    virtual void add_block(const std::unique_ptr<DataType[]>& dat);
    virtual void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0);
    virtual void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1);
    virtual void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2);
    virtual void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3);
    virtual void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                   const Index& i4);
    virtual void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                   const Index& i4, const Index& i5);
    virtual void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                   const Index& i4, const Index& i5, const Index& i6);
    virtual void add_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                   const Index& i4, const Index& i5, const Index& i6, const Index& i7);

    template<typename ...args>
    size_t blocksize(const args& ...p) const {
      auto a = hashtable_.find(generate_hash_key(p...));
      return a != hashtable_.end() ? (a->second.second - a->second.first + 1) : 0lu;
    }
    template<typename ...args>
    size_t blocksize_alloc(const args& ...p) const {
      auto a = hashtable_.find(generate_hash_key(p...));
      return (a != hashtable_.end() && initialized_) ? (a->second.second - a->second.first + 1) : 0lu;
    }

    size_t size() const { return totalsize_; }
    size_t size_alloc() const { return initialized_ ? size() : 0lu; }

    void initialize();

    void zero();
    void scale(const DataType& a);

    template<typename ...args>
    bool is_local(args&& ...p) const { return is_local(generate_hash_key(p...)); }
    bool is_local(const size_t key) const;

    StorageIncore<DataType>& operator=(const StorageIncore<DataType>& o);
    void ax_plus_y(const DataType& a, const StorageIncore<DataType>& o);
    void ax_plus_y(const DataType& a, const std::shared_ptr<StorageIncore<DataType>> o) { ax_plus_y(a, *o); };
    DataType dot_product(const StorageIncore<DataType>& o) const;

    virtual void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& p) { }
};

template<> double StorageIncore<double>::dot_product(const StorageIncore<double>& o) const;
template<> std::complex<double> StorageIncore<std::complex<double>>::dot_product(const StorageIncore<std::complex<double>>& o) const;


extern template class StorageIncore<double>;
extern template class StorageIncore<std::complex<double>>;

template<typename DataType>
using Storage = StorageIncore<DataType>;

}
}

#endif

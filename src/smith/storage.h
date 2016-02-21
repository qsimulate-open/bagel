//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: storage.h
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
#include <list>
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


namespace {
  void arg_convert_impl(std::vector<Index>& a) { }
  template<typename... args>
  void arg_convert_impl(std::vector<Index>& a, const Index& i, args... tail) {
    a.push_back(i);
    arg_convert_impl(a, tail...);
  }
  template<typename... args>
  std::vector<Index> arg_convert(args... p) {
    std::vector<Index> a;
    arg_convert_impl(a, p...);
    return a;
  }
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
    void put_block_(const std::unique_ptr<DataType[]>& dat, const size_t& key);
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

    virtual void put_block(const std::unique_ptr<DataType[]>& dat);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                   const Index& i4);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                   const Index& i4, const Index& i5);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                   const Index& i4, const Index& i5, const Index& i6);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                   const Index& i4, const Index& i5, const Index& i6, const Index& i7);
    virtual void put_block(const std::unique_ptr<DataType[]>& dat, const std::vector<Index> i);

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

    size_t blocksize() const { return 1lu; }
    template<typename ...args>
    size_t blocksize(const Index& i, args&& ...p) const { return i.size()*blocksize(p...); }
    size_t blocksize(const std::vector<Index>& p) const { return std::accumulate(p.begin(), p.end(), 1lu, [](size_t a, const Index& b){ return a*b.size(); }); }

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

    // for Kramers storage
    virtual void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& p) { }
    virtual void set_stored_sectors(const std::list<std::vector<bool>>& p) { }
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

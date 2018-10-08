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

#include <map>
#include <unordered_map>
#include <list>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <src/smith/indexrange.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/parallel/rmawindow.h>

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
class StorageIncore : public RMAWindow<DataType> {
  public:
    using RMAWindow<DataType>::initialize;
    using RMAWindow<DataType>::initialized;
    using RMAWindow<DataType>::rma_get;
    using RMAWindow<DataType>::rma_put;
    using RMAWindow<DataType>::rma_add;

  protected:
    size_t totalsize_;

    // this relates hash keys to lo and high of the block
    std::unordered_map<size_t, std::pair<size_t, size_t>> hashtable_;
    // distribution information. Relates lo and the process number
    std::map<size_t, int> blocks_;

    // local storage
    size_t local_lo_;
    size_t local_hi_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      std::map<size_t, std::pair<size_t, size_t>> hashtable_ordered;
      for (auto& i: hashtable_)
        hashtable_ordered.emplace(i);
      ar << RMAWindow<DataType>::initialized_ << totalsize_ << hashtable_ordered;

      // Process 0 collects and saves tensor's contents, tile by tile
      // Other processes do nothing; this requires them to be writing to a different file from Process 0
      if (mpi__->rank() == 0) {
        for (auto& i : hashtable_ordered) {
          size_t rank, off, size;
          std::tie(rank, off, size) = locate(i.first);
          std::vector<DataType> tmp(size, 0.0);
          rma_get(tmp.data(), i.first);
          ar << tmp;
        }
      }
      mpi__->barrier();
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      std::map<size_t, std::pair<size_t, size_t>> hashtable_ordered;
      bool init;
      ar >> init >> totalsize_ >> hashtable_ordered;

      // Determine distribution information (assuming mpi__->size() might have changed)
      const size_t blocksize = (totalsize_-1)/mpi__->size()+1;
      blocks_ = std::map<size_t, int>();
      size_t tsize = 0;
      for (auto& i : hashtable_ordered) {
        if (i.second.first == i.second.second) continue;
        if (blocks_.size()*blocksize <= tsize)
          blocks_.emplace(tsize, blocks_.size());
        tsize += (i.second.second - i.second.first);
      }
      assert(totalsize_ == tsize);
      local_lo_ = local_hi_ = std::numeric_limits<size_t>::max();
      for (auto iter = blocks_.begin(); iter != blocks_.end(); ++iter)
        if (iter->second == mpi__->rank()) {
          local_lo_ = iter->first;
          auto iterp = ++iter;
          local_hi_ = iterp != blocks_.end() ? iterp->first : totalsize_;
          break;
        }

      // copy hashtable
      for (auto& i: hashtable_ordered)
        hashtable_.emplace(i);
      if (init)
        initialize();

      // All processes read the whole archive, and save the data that belong to them
      for (auto& i : hashtable_ordered) {
        size_t rank, off, size;
        std::tie(rank, off, size) = locate(i.first);
        std::vector<DataType> tmp(size, 0.0);
        ar >> tmp;
        if (rank == mpi__->rank())
          rma_put(tmp.data(), i.first);
      }
      mpi__->barrier();
    }

  public:
    StorageIncore() { }
    StorageIncore(const std::map<size_t, size_t>& size, bool init);

    // required functions by RMAWindow
    bool is_local(const size_t key) const override;
    size_t localsize() const override;
    std::tuple<size_t, size_t, size_t> locate(const size_t key) const override; // returns (process, offset, size)

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
    size_t size_alloc() const { return initialized() ? size() : 0lu; }

    template<typename ...args>
    bool is_local(args&& ...p) const { return is_local(generate_hash_key(p...)); }

    // for Kramers storage
    virtual void set_perm(const std::map<std::vector<int>, std::pair<double,bool>>& p) { }
    virtual void set_stored_sectors(const std::list<std::vector<bool>>& p) { }
};

extern template class StorageIncore<double>;
extern template class StorageIncore<std::complex<double>>;

template<typename DataType>
using Storage = StorageIncore<DataType>;

}
}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SMITH::StorageIncore<double>)
BOOST_CLASS_EXPORT_KEY(bagel::SMITH::StorageIncore<std::complex<double>>)

#endif

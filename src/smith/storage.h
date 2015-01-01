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
#include <memory>
#include <tuple>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <fstream>
#include <cstdio>

namespace bagel {
namespace SMITH {

class StorageBlock {
  protected:
    std::unique_ptr<double[]> data_;
    size_t size_;
    bool initialized_;

    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }
  public:
    StorageBlock(const size_t size, const bool init) : size_(size), initialized_(init) {
      if (init) {
        data_ = std::unique_ptr<double[]>(new double[size_]);
        zero();
      }
    }

    void zero() {
      if (initialized_)
        std::fill_n(data(), size_, 0.0);
    }

    size_t size() const { return size_; }
    size_t size_alloc() const { return initialized_ ? size_ : 0lu; }

    StorageBlock& operator=(const StorageBlock& o) {
      if (o.initialized_ && !initialized_) {
        data_ = std::unique_ptr<double[]>(new double[size_]);
        initialized_ = true;
      }
      if (o.initialized_)
        std::copy_n(o.data(), size_, data());
      return *this;
    }

    void put_block(std::unique_ptr<double[]>&& o) {
      assert(!initialized_);
      initialized_ = true;
      data_ = std::move(o);
    }

    std::unique_ptr<double[]> get_block() const {
      assert(initialized_);
      std::unique_ptr<double[]> out(new double[size_]);
      std::copy_n(data_.get(), size_, out.get());
      return std::move(out);
    }

    std::unique_ptr<double[]> move_block() {
      if (!initialized_) {
        initialized_ = true;
        data_ = std::unique_ptr<double[]>(new double[size_]);
        zero();
      }
      initialized_ = false;
      return std::move(data_);
    }

    void add_block(const std::unique_ptr<double[]>& o) {
      assert(initialized_);
      blas::ax_plus_y_n(1.0, o.get(), size_, data());
    }

    double dot_product(const StorageBlock& o) const {
      assert(size_ == o.size_ && !(initialized_ ^ o.initialized_));
      return initialized_ ? blas::dot_product(data(), size_, o.data()) : 0.0;
    }

    void ax_plus_y(const double a, const StorageBlock& o) {
      assert(size_ == o.size_ && !(initialized_ ^ o.initialized_));
      if (initialized_)
        blas::ax_plus_y_n(a, o.data(), size_, data());
    }

    void scale(const double a) {
      if (initialized_)
        blas::scale_n(a, data(), size_);
    }
};

class Storage_base {
  protected:
    // this relates hash keys, block number, and block lengths (in this order).
    std::map<size_t, std::shared_ptr<StorageBlock>> hashtable_;

  public:
    // size contains hashkey and length (in this order)
    Storage_base(const std::map<size_t, size_t>& size, bool init) {
      for (auto& i : size)
        hashtable_.emplace(i.first, std::make_shared<StorageBlock>(i.second, init));
    }

    // functions that return protected members
    size_t blocksize(const size_t hash) const {
      auto a = hashtable_.find(hash);
      return a != hashtable_.end() ? a->second->size() : 0lu;
    }
    size_t blocksize_alloc(const size_t hash) const {
      auto a = hashtable_.find(hash);
      return a != hashtable_.end() ? a->second->size_alloc() : 0lu;
    }


    // get, move, put, and add a block from the storage and returns unique_ptr<double[]>, which is local
    virtual std::unique_ptr<double[]> get_block(const size_t& key) const = 0;
    virtual std::unique_ptr<double[]> move_block(const size_t& key) = 0;
    virtual void put_block(const size_t& key, std::unique_ptr<double[]>& dat) = 0;
    virtual void add_block(const size_t& key, const std::unique_ptr<double[]>& dat) = 0;

    virtual void zero() = 0;
    virtual void scale(const double a) = 0;

};

class Storage_Incore : public Storage_base {
  public:
    Storage_Incore(const std::map<size_t, size_t>& size, bool init);

    std::unique_ptr<double[]> get_block(const size_t& key) const;
    std::unique_ptr<double[]> move_block(const size_t& key);
    void put_block(const size_t& key, std::unique_ptr<double[]>& dat);
    void add_block(const size_t& key, const std::unique_ptr<double[]>& dat);

    void zero();
    void scale(const double a);

    Storage_Incore& operator=(const Storage_Incore& o);
    void ax_plus_y(const double a, const Storage_Incore& o);
    void ax_plus_y(const double a, const std::shared_ptr<Storage_Incore> o) { ax_plus_y(a, *o); };
    double dot_product(const Storage_Incore& o) const;

};

#if 0
class Storage_Disk : public Storage_base {
  protected:
    std::string filename_;
    mutable std::fstream data_;

    size_t totalsize_;
    std::map<size_t, size_t> offset_;

    const long long cachesize_ = 100000lu;

  public:
    Storage_Disk(const std::map<size_t, size_t>& size, bool init);
    ~Storage_Disk() { std::remove(filename_.c_str()); }

    std::unique_ptr<double[]> get_block(const size_t& key) const;
    std::unique_ptr<double[]> move_block(const size_t& key);
    void put_block(const size_t& key, std::unique_ptr<double[]>& dat);
    void add_block(const size_t& key, const std::unique_ptr<double[]>& dat);

    void zero();
    void scale(const double a);

    Storage_Disk& operator=(const Storage_Disk& o);
    void ax_plus_y(const double a, const Storage_Disk& o);
    void ax_plus_y(const double a, const std::shared_ptr<Storage_Disk> o) { ax_plus_y(a, *o); };
    double dot_product(const Storage_Disk& o) const;

};
#endif

//#ifdef SMITH_INCORE
#if 1
using Storage = Storage_Incore;
#else
using Storage = Storage_Disk;
#endif

}
}

#endif

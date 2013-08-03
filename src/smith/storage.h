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
// the Free Software Foundation; either version 2, or (at your option)
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

namespace bagel {
namespace SMITH {

class Storage_base {
  protected:
    // length of this storage
    size_t length_;
    // this relates hash keys, block number, and block lengths (in this order).
    std::map<size_t, std::pair<size_t, size_t>> hashtable_;
    std::vector<bool> initialized_;


  public:
    // size contains hashkey and length (in this order)
    Storage_base(const std::map<size_t, size_t>& size, bool init) {
      length_ = 0lu;
      size_t cnt = 0;
      for (auto i = size.begin(); i != size.end(); ++i, ++cnt) {
        auto j = hashtable_.insert(std::make_pair(i->first, std::make_pair(cnt, i->second)));
        if (!j.second) throw std::logic_error("duplicated hash keys in Storage::Storage");
        length_ += i->second;
      }
      initialized_ = std::vector<bool>(cnt, init);
    }

    // functions that return protected members
    size_t length() const { return length_; }
    size_t blocksize(const size_t hash) const {
      auto a = hashtable_.find(hash);
      return a != hashtable_.end() ? a->second.second : 0;
    }

    // get, move, put, and add a block from the storage and returns unique_ptr<double[]>, which is local
    virtual std::unique_ptr<double[]> get_block(const size_t& key) const = 0;
    virtual std::unique_ptr<double[]> move_block(const size_t& key) = 0;
    virtual void put_block(const size_t& key, std::unique_ptr<double[]>& dat) = 0;
    virtual void add_block(const size_t& key, const std::unique_ptr<double[]>& dat) = 0;

    virtual void zero() = 0;
    virtual void scale(const double a) = 0;

    virtual void initialize() = 0;

    bool initialized(const int i) const { return initialized_[i]; }

};

class Storage_Incore : public Storage_base {
  protected:
    std::vector<std::unique_ptr<double[]>> data_;

  public:
    Storage_Incore(const std::map<size_t, size_t>& size, bool init);

    std::unique_ptr<double[]> get_block(const size_t& key) const;
    std::unique_ptr<double[]> move_block(const size_t& key);
    void put_block(const size_t& key, std::unique_ptr<double[]>& dat);
    void add_block(const size_t& key, const std::unique_ptr<double[]>& dat);

    void zero();
    void scale(const double a);

    Storage_Incore& operator=(const Storage_Incore& o);
    void daxpy(const double a, const Storage_Incore& o);
    void daxpy(const double a, const std::shared_ptr<Storage_Incore> o) { daxpy(a, *o); };
    double ddot(const Storage_Incore& o) const;

    void initialize();
};

}
}

#endif

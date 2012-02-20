//
// Newint - Parallel electron correlation program.
// Filename: indexrange.h
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


#ifndef __SRC_SMITH_INDEXRANGE_H
#define __SRC_SMITH_INDEXRANGE_H

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace SMITH {

// one index block
class Index {
  protected:
    size_t offset_;
    size_t size_;
    size_t key_;
  public:
    Index(const size_t& o, const size_t& s, const size_t& i) : offset_(o), size_(s), key_(i) {};
    ~Index() {};
    size_t offset() const { return offset_; };
    size_t size() const { return size_; };
    size_t key() const { return key_; };

    bool operator==(const Index& o) {
      return offset_ == o.offset_ && size_ == o.size_ && key_ == o.key_;
    };
};

class IndexRange {
  protected:
    // pair of offset and size (in this order)
    std::vector<Index> range_;

    // total size of this index range
    int size_;
    int keyoffset_;

  public:
    IndexRange(const int size, const int maxblock = 10, const int boffset = 0) : keyoffset_(boffset) {
      // first determine number of blocks. 
      const size_t nbl = (size-1) / maxblock + 1;
      const size_t nblock = (size-1) / nbl + 1;
      // we want to distribute orbitals as evenly as possible 
      const size_t rem = nbl * nblock - size; 
      std::vector<size_t> blocksizes(nbl, nblock); 
      auto iter = blocksizes.rbegin();
      for (int k = 0; k != rem; ++iter, ++k) --*iter;
      // push back to range_
      size_t off = 0;
      // key is offsetted by the input value
      size_t cnt = boffset;
      for (auto i = blocksizes.begin(); i != blocksizes.end(); ++i, ++cnt) {
        Index t(off, *i, cnt);
        range_.push_back(t);
        off += *i;
      }
      // set size_
      size_ = off;
    }; 
    ~IndexRange() {};

    const std::vector<Index>& range() const { return range_; };
    Index range(const int i) const { return range_[i]; };

    int nblock() const { return range_.size(); };
    int size() const { return size_; };
    int keyoffset() const { return keyoffset_; };

    bool operator==(const IndexRange& o) {
      bool out = size_ == o.size_;
      if (range_.size() == o.range_.size()) {
        auto i = range_.begin();
        for (auto j = o.range_.begin(); i != range_.end(); ++i, ++j) out &= (*i) == (*j);
      } else {
        out = false;
      }
      return out;
    };
    bool operator!=(const IndexRange& o) { return !(*this == o); };

    std::string str() const {
      std::stringstream ss;
      for (auto i = range_.begin(); i != range_.end(); ++i)
        ss << std::setw(10) << i->offset() << std::setw(10) << i->size() << std::endl; 
      return ss.str();
    };
    void print() const { std::cout << str() << std::endl; };
};

}

#endif

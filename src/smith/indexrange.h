//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: indexrange.h
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


#ifndef __SRC_SMITH_INDEXRANGE_H
#define __SRC_SMITH_INDEXRANGE_H

#include <stddef.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <src/util/parallel/staticdist.h>
#include <src/util/serialization.h>

namespace bagel {
namespace SMITH {

// one index block
class Index {
  protected:
    size_t offset_;
    size_t offset2_;
    size_t size_;
    size_t key_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & offset_ & offset2_ & size_ & key_;
    }

  public:
    Index(const size_t& o, const size_t& o2, const size_t& s, const size_t& i) : offset_(o), offset2_(o2), size_(s), key_(i) {}
    Index() {}
    ~Index() {}
    size_t offset() const { return offset_; }
    size_t size() const { return size_; }
    size_t key() const { return key_; }

    // returns if this block corresponds to Kramers + or -.
    bool kramers() const { return offset_ == offset2_; }
    // returns the offset of corresponding Kramers+ Index
    size_t kramers_offset() const { return offset2_; }

    bool operator==(const Index& o) const {
      return offset_ == o.offset_ && offset2_ == o.offset2_ && size_ == o.size_ && key_ == o.key_;
    }
};

class IndexRange {
  protected:
    // pair of offset and size (in this order)
    std::vector<Index> range_;

    // total size of this index range
    int size_;
    int keyoffset_;
    int orboffset_;

    // set to be an offset for the spin/Kramers counterpart if necessary
    int orboffset2_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & range_ & size_ & keyoffset_ & orboffset_ & orboffset2_;
    }

  public:
    IndexRange(const int size, const int maxblock = 10, const int boffset = 0, const int orboff = 0, const int orboff2 = -1);

    // make IndexRange based on StaticDist. The offset in StaticDist is ignored
    IndexRange(std::shared_ptr<const StaticDist> dist, const int boffset = 0, const int orboff = 0, const int orboff2 = -1);

    IndexRange() {}
    ~IndexRange() {}

    const std::vector<Index>& range() const { return range_; }
    Index range(const int i) const { return range_[i]; }

    std::vector<Index>::iterator begin() { return range_.begin(); }
    std::vector<Index>::iterator end() { return range_.end(); }
    std::vector<Index>::const_iterator begin() const { return range_.begin(); }
    std::vector<Index>::const_iterator end() const { return range_.end(); }
    Index& front() { return range_.front(); }
    Index& back() { return range_.back(); }
    const Index& front() const { return range_.front(); }
    const Index& back() const { return range_.back(); }

    int nblock() const { return range_.size(); }
    int size() const { return size_; }
    int keyoffset() const { return keyoffset_; }

    void merge(const IndexRange& o);

    bool operator==(const IndexRange& o) const;
    bool operator!=(const IndexRange& o) const { return !(*this == o); }

    std::string str() const;
    void print() const { std::cout << str() << std::endl; }
};

}
}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SMITH::Index)
BOOST_CLASS_EXPORT_KEY(bagel::SMITH::IndexRange)

#endif

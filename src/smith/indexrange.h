//
// Author : Toru Shiozaki
// Date   : Feb 2012
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
};

class IndexRange {
  protected:
    // pair of offset and size (in this order)
    std::vector<Index> range_;

    // total size of this index range
    int size_;

  public:
    IndexRange(const int size, const int maxblock = 10) {
      // first determine number of blocks. 
      const size_t nbl = (size-1) / maxblock + 1;
      const size_t nblock = size / nbl + 1;
      // we want to distribute orbitals as evenly as possible 
      const size_t rem = nbl * nblock - size; 
      std::vector<size_t> blocksizes(nbl, nblock); 
      auto iter = blocksizes.rbegin();
      for (int k = 0; k != rem; ++iter, ++k) --*iter;
      // push back to range_
      size_t off = 0;
      size_t cnt = 0;
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

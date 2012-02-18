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

class IndexRange {
  protected:
    // pair of offset and size (in this order)
    std::list<std::pair<size_t, size_t> > range_;

    // total size of this index range
    int size_;

  public:
    IndexRange(const int size, const int maxblock = 10) {
      // first determine number of blocks. 
      const size_t nbl = (size-1) / maxblock + 1;
      // we want to distribute orbitals as evenly as possible 
      const size_t rem = nbl *maxblock - size; 
      std::vector<size_t> blocksizes(nbl, maxblock); 
      auto iter = blocksizes.rbegin();
      for (int k = 0; k != rem; ++iter, ++k) --*iter;
      // push back to range_
      int off = 0;
      for (auto i = blocksizes.begin(); i != blocksizes.end(); ++i) {
        range_.push_back(std::make_pair(off, *i));
        off += *i;
      }
      // set size_
      size_ = off;
    }; 
    ~IndexRange() {};

    const std::list<std::pair<size_t, size_t> >& range() const { return range_; };
    int nblock() const { return range_.size(); };
    int size() const { return size_; };

    std::string str() const {
      std::stringstream ss;
      for (auto i = range_.begin(); i != range_.end(); ++i)
        ss << std::setw(10) << i->first << std::setw(10) << i->second << std::endl; 
      return ss.str();
    };
    void print() const { std::cout << str() << std::endl; };
};

}

#endif

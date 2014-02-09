//
// BAGEL - Parallel electron correlation program.
// Filename: cistringset.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@u.northwestern.edu>
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

#ifndef __SRC_CIUTIL_CISTRINGSET_H
#define __SRC_CIUTIL_CISTRINGSET_H

#include <src/ciutil/cistring.h>

namespace bagel {

// A set of CIStrings (implemented as a list) that has the same #electrons
// which appears often in the RAS implementation
// Interface is assumed to be identical to that of CIStrings

template <class StringType>
class CIStringSet {
  protected:
    std::map<size_t, std::shared_ptr<StringType>> stringset_;

    // global information
    int nele_;
    int norb_;
    size_t size_;

    std::vector<std::bitset<nbit__>> strings_;

    std::shared_ptr<StringMap> phi_;

  public:
    CIStringSet(const std::list<std::shared_ptr<StringType>>& o) {
      // copy construct with an offset
      nele_ = o.front()->nele();
      norb_ = o.front()->norb();
      size_ = 0ull;
      for (auto& i : o) {
        // checks if CIString has the same #electrons and #orbitals
        assert(nele_ == i->nele() && norb_ == i->norb());
        auto tmp = std::make_shared<StringType>(*i, size_);
        stringset_.push_back(tmp);
        size_ += tmp->size();
        strings_.insert(strings_.end(), tmp->strings().begin(), tmp->strings().end());
      }
      assert(size_ == strings_.size());
    }

    int nele() const { return nele_; }
    int norb() const { return norb_; }
    size_t size() const { return size_; }

    const std::vector<std::bitset<nbit__>>& strings() const { return strings_; }
    const std::bitset<nbit__>& strings(const size_t i) const { return strings_[i]; }

    std::vector<std::bitset<nbit__>>::iterator begin() { return strings_.begin(); }
    std::vector<std::bitset<nbit__>>::iterator end() { return strings_.end(); }
    std::vector<std::bitset<nbit__>>::const_iterator begin() const { return strings_.cbegin(); }
    std::vector<std::bitset<nbit__>>::const_iterator end() const { return strings_.cend(); }

};

}

#endif

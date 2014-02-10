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
#include <src/ciutil/bitutil.h>

namespace bagel {

// A set of CIStrings (implemented as a list) that has the same #electrons
// which appears often in the RAS implementation
// Interface is assumed to be identical to that of CIStrings

template <class StringType>
class CIStringSet {
  protected:
    std::list<std::shared_ptr<StringType>> stringset_;

    // global information
    int nele_;
    int norb_;
    size_t size_;

    std::vector<std::bitset<nbit__>> strings_;

    std::shared_ptr<StringMap> phi_;

    void init() {
      construct_phi();
    }
    void construct_phi() {
      phi_ = std::make_shared<StringMap>(size_);
      phi_->reserve(norb_*norb_);

      std::unordered_map<size_t, size_t> lexmap;
      for (size_t i = 0; i < size_; ++i)
        lexmap[strings_[i].to_ullong()] = i;

      size_t tindex = 0;
      for (auto& istring : strings_) {
        for (int j = 0; j < norb_; ++j) {
          if (!istring[j]) continue;
          std::bitset<nbit__> intermediatebit = istring; intermediatebit.reset(j);
          for (int i = 0; i < norb_; ++i) {
            if (intermediatebit[i]) continue;
            std::bitset<nbit__> sourcebit = intermediatebit; sourcebit.set(i);
            if (allowed(sourcebit)) {
              assert(lexmap.find(sourcebit.to_ullong()) != lexmap.end());
              (*phi_)[tindex].emplace_back(tindex, sign(istring, i, j), lexmap[sourcebit.to_ullong()], j+i*norb_);
            }
          }
        }
        (*phi_)[tindex++].shrink_to_fit();
      }
    }

  public:
    CIStringSet(const std::list<std::shared_ptr<const StringType>>& o) {
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

      init();
    }

    int nele() const { return nele_; }
    int norb() const { return norb_; }
    size_t size() const { return size_; }
    size_t nspaces() const { return stringset_.size(); }

    bool allowed(const std::bitset<nbit__>& bit) const {
      return std::any_of(stringset_.begin(), stringset_.end(),
                         [&](const std::shared_ptr<StringType>& a) { return a->contains(bit); });
    }

    const std::vector<std::bitset<nbit__>>& strings() const { return strings_; }
    const std::bitset<nbit__>& strings(const size_t i) const { return strings_[i]; }

    typename std::list<std::shared_ptr<StringType>>::iterator begin() { return stringset_.begin(); }
    typename std::list<std::shared_ptr<StringType>>::iterator end() { return stringset_.end(); }
    typename std::list<std::shared_ptr<StringType>>::const_iterator begin() const { return stringset_.cbegin(); }
    typename std::list<std::shared_ptr<StringType>>::const_iterator end() const { return stringset_.cend(); }

    std::shared_ptr<const StringMap> phi() const { return phi_; }

    std::shared_ptr<const StringType> find_string(const std::bitset<nbit__>& bit) const {
      for (auto& i : stringset_)
        if (i->contains(bit)) return i;
      return  std::make_shared<const StringType>();
    }

    std::shared_ptr<const StringType> find_string(const int nholes, const int nparticles) const {
      for (auto& i : stringset_)
        if (i->matches(nholes, nparticles)) return i;
      return  std::make_shared<const StringType>();
    }

    size_t lexical_offset(const std::bitset<nbit__>& bit) const {
      return find_string(bit)->lexical_offset(bit);
    }

    size_t lexical_zero(const std::bitset<nbit__>& bit) const {
      return find_string(bit)->lexical_zero(bit);
    }
};

}

#endif

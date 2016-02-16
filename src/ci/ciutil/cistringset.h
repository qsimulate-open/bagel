//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cistringset.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@u.northwestern.edu>
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

#ifndef __SRC_CIUTIL_CISTRINGSET_H
#define __SRC_CIUTIL_CISTRINGSET_H

#include <src/ci/ciutil/citraits.h>
#include <src/ci/ciutil/cistring.h>
#include <src/ci/ciutil/bitutil.h>

namespace bagel {

// A set of CIStrings (implemented as a list) that has the same #electrons
// which appears often in the RAS implementation
// Interface is assumed to be identical to that of CIStrings

template <class StringType,
          class = typename std::enable_if<is_cistring<StringType>::value>::type
         >
class CIStringSet {
  protected:
    std::list<std::shared_ptr<StringType>> stringset_;

    // global information
    int nele_;
    int norb_;
    size_t size_;

    std::vector<std::bitset<nbit__>> strings_;

    std::shared_ptr<StringMap> phi_;
    std::shared_ptr<StringMap> uncompressed_phi_;

    void init() {
      construct_phi();
    }
    void construct_phi(); // should be specialized

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & stringset_ & nele_ & norb_ & size_ & strings_ & phi_ & uncompressed_phi_;
    }

  public:
    CIStringSet() { }
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

    // use nbit__ for the striding since nbit__ should be an upper bound on nele_ anyways
    size_t key() const { return nele_+nbit__*norb_; }

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
    std::shared_ptr<const StringMap> uncompressed_phi() const { return uncompressed_phi_; }

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

template<>
void CIStringSet<FCIString>::construct_phi();
template<>
void CIStringSet<RASString>::construct_phi();

}

extern template class bagel::CIStringSet<bagel::FCIString>;
extern template class bagel::CIStringSet<bagel::RASString>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::CIStringSet<bagel::FCIString>)
BOOST_CLASS_EXPORT_KEY(bagel::CIStringSet<bagel::RASString>)

#endif

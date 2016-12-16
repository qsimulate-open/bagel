//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cistringspace.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_CIUTIL_CISTRINGSPACE_H
#define __SRC_CIUTIL_CISTRINGSPACE_H

#include <unordered_map>
#include <src/util/serialization.h>
#include <src/ci/ciutil/cistringset.h>

namespace bagel {

template <class StringSetType,
          class = typename std::enable_if<is_cistringset<StringSetType>::value>::type
         >
class CIStringSpace {
  protected:
    std::unordered_map<int, std::shared_ptr<const StringSetType>> strings_;

    // it is assumed that every StringSetType has the same norb_.
    size_t norb_;

    std::unordered_map<size_t, std::shared_ptr<StringMap>> phiup_;
    std::unordered_map<size_t, std::shared_ptr<StringMap>> phidown_;

    static size_t hash(const std::shared_ptr<const StringSetType>& a) { return a->key(); }

    void build_linkage(std::shared_ptr<const StringSetType> a, std::shared_ptr<const StringSetType> b, const int fac) {
      std::shared_ptr<const StringSetType> plus, ref;

      const int diff = a->nele() - b->nele();
      if (diff == 1)
        std::tie(ref, plus) = std::make_pair(b, a);
      else if (diff == -1)
        std::tie(ref, plus) = std::make_pair(a, b);
      else
        // quick return if a and b cannot be connected by a single excitation.
        return;

      auto phiup   = std::make_shared<StringMap>(norb_);
      auto phidown = std::make_shared<StringMap>(norb_);
      phiup->reserve(plus->size());
      phidown->reserve(ref->size());

      std::vector<std::bitset<nbit__>> stringplus = plus->strings();
      std::vector<std::bitset<nbit__>> string     = ref->strings();


      for (auto& istring : string) {
        for (size_t i = 0; i != norb_; ++i) {
          if (!istring[i]) { // creation
            const size_t source = ref->lexical_offset(istring);
            std::bitset<nbit__> nbit = istring; nbit.set(i); // created.
            if (plus->allowed(nbit)) {
              const size_t target = plus->lexical_offset(nbit);
              const int s = sign(nbit, i) * fac;
              (*phiup)[i].emplace_back(target, s, source, 0);
              (*phidown)[i].emplace_back(source, s, target, 0);
            }
          }
        }
      }

      phiup->shrink_to_fit();
      phidown->shrink_to_fit();

      const size_t upkey = ref->key();
      if (phiup_.find(upkey) == phiup_.end())
        phiup_[upkey] = phiup;
      else
        phiup_[upkey]->insert(phiup);

      const size_t downkey = plus->key();
      if (phidown_.find(downkey) == phidown_.end())
        phidown_[downkey] = phidown;
      else
        phidown_[downkey]->insert(phidown);
    }

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & strings_ & norb_ & phiup_ & phidown_;
    }

  public:
    CIStringSpace() { }

    CIStringSpace(std::initializer_list<std::shared_ptr<const StringSetType>> s)
      : CIStringSpace(std::list<std::shared_ptr<const StringSetType>>(s.begin(), s.end())) {
    }

    CIStringSpace(std::list<std::shared_ptr<const StringSetType>> s) {
      assert(!s.empty());
      norb_ = s.front()->norb();
      for (auto& i : s) {
        if (norb_ != i->norb())
          throw std::logic_error("All CIStrings in CIStringSpace should have the same norb.");
        strings_.emplace(i->nele(), i);
      }
    }

    void build_linkage(const int fac = 1) {
      // every time this function is called, we construct phiup and down from scratch (cheap anyway)
      phiup_.clear();
      phidown_.clear();
      for (auto it = strings_.begin(); it != strings_.end(); ++it) {
        auto jt = it; ++jt;
        for ( ; jt != strings_.end(); ++jt)
          build_linkage(it->second, jt->second, fac);
      }
    }

    std::shared_ptr<const StringMap> phiup(const std::shared_ptr<const StringSetType>& a) const {
      return phiup_.find(hash(a))->second;
    }

    std::shared_ptr<const StringMap> phidown(const std::shared_ptr<const StringSetType>& a) const {
      return phidown_.find(hash(a))->second;
    }

};

}

#endif

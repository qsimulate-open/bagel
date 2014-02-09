//
// BAGEL - Parallel electron correlation program.
// Filename: cistringspace.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_CIUTIL_CISTRINGSPACE_H
#define __SRC_CIUTIL_CISTRINGSPACE_H

#include <unordered_map>
#include <src/util/serialization.h>
#include <src/ciutil/cistringmap.h>
#include <src/ciutil/bitutil.h>
#include <src/ciutil/cistring.h>

namespace bagel {

template <class StringType>
class CIStringSpace {
  protected:
    std::list<std::shared_ptr<const StringType>> strings_;

    // it is assumed that every StringType has the same norb_.
    size_t norb_;

    std::unordered_map<size_t, std::shared_ptr<StringMap>> phiup_;
    std::unordered_map<size_t, std::shared_ptr<StringMap>> phidown_;

    static size_t hash(const std::shared_ptr<const StringType>& a) { return a->key(); }

    void build_linkage(std::shared_ptr<const StringType> a, std::shared_ptr<const StringType> b, const int fac) {
      std::shared_ptr<const StringType> plus, ref;

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
            if (plus->contains(nbit)) {
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

  public:
    CIStringSpace(std::initializer_list<std::shared_ptr<const StringType>> s)
      : CIStringSpace(std::list<std::shared_ptr<const StringType>>(s.begin(), s.end())) {
    }

    CIStringSpace(std::list<std::shared_ptr<const StringType>> s) : strings_(s) {
      assert(!strings_.empty());
      norb_ = strings_.front()->norb();
      for (auto& i : strings_)
        if (norb_ != i->norb())
          throw std::logic_error("All CIStrings in CIStringSpace should have the same norb.");
    }

    void build_linkage(const int fac = 1) {
      // every time this function is called, we construct phiup and down from scratch (cheap anyway)
      phiup_.clear();
      phidown_.clear();
      for (auto it = strings_.begin(); it != strings_.end(); ++it) {
        auto jt = it; ++jt;
        for ( ; jt != strings_.end(); ++jt)
          build_linkage(*it, *jt, fac);
      }
    }

    std::shared_ptr<const StringMap> phiup(const std::shared_ptr<const StringType>& a) const {
      return phiup_.find(hash(a))->second;
    }

    std::shared_ptr<const StringMap> phidown(const std::shared_ptr<const StringType>& a) const {
      return phidown_.find(hash(a))->second;
    }

};

}

#endif

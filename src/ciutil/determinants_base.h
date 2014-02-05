//
// BAGEL - Parallel electron correlation program.
// Filename: determinants_base.h
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


#ifndef __SRC_CIUTIL_DETERMINANTS_BASE_H
#define __SRC_CIUTIL_DETERMINANTS_BASE_H

#include <tuple>
#include <string>
#include <map>
#include <algorithm>
#include <cassert>
#include <src/util/constants.h>
#include <src/util/serialization.h>
#include <src/ciutil/cistring.h>

namespace bagel {

// Determinants class without linkages.
class Determinants_base {
  protected:
    // assuming that the number of active orbitals are the same in alpha and beta.
    int norb_;

    int nelea_;
    int neleb_;

    // string lists
    std::shared_ptr<FCIString> astring_;
    std::shared_ptr<FCIString> bstring_;

    // this is slow but robust implementation of bit to number converter.
    std::vector<int> bit_to_numbers(std::bitset<nbit__> bit) const {
      std::vector<int> out;
      for (int i = 0; i != norb_; ++i) if (bit[i]) out.push_back(i);
      return out;
    }

    std::bitset<nbit__> numbers_to_bit(const std::vector<int>& num) const {
      std::bitset<nbit__> out(0);
      for (auto& i : num) out.set(i);
      return out;
    }

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & norb_ & nelea_ & neleb_ & astring_ & bstring_;
    }

  public:
    Determinants_base() { }
    Determinants_base(const int norb, const int nelea, const int neleb, const bool mute = false);
    virtual ~Determinants_base() { }

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    // string size
    std::tuple<size_t, size_t> len_string() const { return std::make_tuple(lena(), lenb()); }

    size_t lena() const { return astring_->size(); }
    size_t lenb() const { return bstring_->size(); }

    size_t ncsfs() const;

    std::string print_bit(std::bitset<nbit__> bit) const {
      std::string out;
      for (int i = 0; i != norb_; ++i) { out += bit[i] ? "1" : "."; }
      return out;
    }

    std::string print_bit(std::bitset<nbit__> bit1, std::bitset<nbit__> bit2) const {
      std::string out;
      for (int i = 0; i != norb_; ++i) {
        if (bit1[i] && bit2[i]) { out += "2"; }
        else if (bit1[i]) { out += "a"; }
        else if (bit2[i]) { out += "b"; }
        else { out += "."; }
      }
      return out;
    }

    template<int spin>
    int sign(std::bitset<nbit__> bit, int i) const {
      static_assert(nbit__ <= sizeof(unsigned long long)*8, "verify Determinant::sign (and other functions)");
      static_assert(spin == 1 || spin == 0, "illegal call of Determinants_base::sign");
      bit &= (1ull << i) - 1ull;
      return (1 - (((bit.count() + spin*nelea_) & 1 ) << 1));
    }

    static int sign(std::bitset<nbit__> bit, int i, int j) {
      // masking irrelevant bits
      int min, max;
      std::tie(min,max) = std::minmax(i,j);
      bit &= ~((1ull << (min+1)) - 1ull);
      bit &= (1ull << max) - 1ull;
      return 1 - ((bit.count() & 1) << 1);
    }

    // maps bit to lexical numbers.
    template <int spin> unsigned int lexical(std::bitset<nbit__> bit) const {
      return spin == 0 ? astring_->lexical(bit) : bstring_->lexical(bit);
    }

    const std::bitset<nbit__>& stringa(int i) const { return astring_->strings(i); }
    const std::bitset<nbit__>& stringb(int i) const { return bstring_->strings(i); }
    const std::vector<std::bitset<nbit__>>& stringa() const { return astring_->strings(); }
    const std::vector<std::bitset<nbit__>>& stringb() const { return bstring_->strings(); }

    std::pair<std::vector<std::tuple<int, int, int>>, double> spin_adapt(const int, std::bitset<nbit__>, std::bitset<nbit__>) const;

    int nspin() const { return nelea_ - neleb_; }
    int norb()  const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Determinants_base)
namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<Determinants_base, T>::value>::type> {
    typedef Determinants_base type;
  };
}

#endif

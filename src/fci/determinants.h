//
// BAGEL - Parallel electron correlation program.
// Filename: determinants.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_FCI_DETERMINANTS_H
#define __SRC_FCI_DETERMINANTS_H

#include <memory>
#include <tuple>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <bitset>
#include <algorithm>
#include <cassert>

#include <src/util/constants.h>

namespace bagel {

struct DetMap {
  unsigned int target;
  unsigned int source;
  int sign;
  DetMap(unsigned int t, int si, unsigned int s) : target(t), source(s), sign(si) {};
};

// implements a determinant space
class Determinants : public std::enable_shared_from_this<Determinants> {
  protected:
    // assuming that the number of active orbitals are the same in alpha and beta.
    const int norb_;

    const int nelea_;
    const int neleb_;

    const bool compress_;

    /* Links to other determinant spaces accessible by one annihilation or creation operation */
    std::weak_ptr<Determinants> detaddalpha_;
    std::weak_ptr<Determinants> detaddbeta_;
    std::weak_ptr<Determinants> detremalpha_;
    std::weak_ptr<Determinants> detrembeta_;

    // Knowles & Handy lexical mapping
    std::vector<unsigned int> zkl_; // contains zkl (Two dimenional array. See the public function).
    // string lists
    std::vector<std::bitset<nbit__>> stringa_;
    std::vector<std::bitset<nbit__>> stringb_;

    // lexical maps (Zkl)
    void const_lexical_mapping_();
    // alpha and beta string lists
    void const_string_lists_();
    // single displacement vectors Phi's
    template <int>
    void const_phis_(const std::vector<std::bitset<nbit__>>& string,
                     std::vector<std::vector<DetMap>>& target, std::vector<std::vector<DetMap>>& uncompressed_target);

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

    // some utility functions
    unsigned int& zkl(int i, int j, int spin) { return zkl_[i*norb_+j+spin*nelea_*norb_]; }
    const unsigned int& zkl(int i, int j, int spin) const { return zkl_[i*norb_+j+spin*nelea_*norb_]; }

    // configuration list i^dagger j compressed
    std::vector<std::vector<DetMap>> phia_;
    std::vector<std::vector<DetMap>> phib_;

    // configuration list i^dagger j uncompressed
    std::vector<std::vector<DetMap>> phia_uncompressed_;
    std::vector<std::vector<DetMap>> phib_uncompressed_;

    // configuration list i^dagger
    std::vector<std::vector<DetMap>> phiupa_;
    std::vector<std::vector<DetMap>> phiupb_;

    // configuration list i
    std::vector<std::vector<DetMap>> phidowna_;
    std::vector<std::vector<DetMap>> phidownb_;

  public:
    Determinants(const int norb, const int nelea, const int neleb, const bool compress = true, const bool mute = false);
    Determinants(std::shared_ptr<const Determinants> o, const bool compress = true, const bool mute = false) :
      Determinants(o->norb(), o->nelea(), o->neleb(), compress, mute) {} // Shortcut to change compression of Det

    // Shortcut to make an uncompressed and muted Determinants with specified # of electrons (used for compatibility with RASDet)
    std::shared_ptr<Determinants> clone(const int nelea, const int neleb) const { return std::make_shared<Determinants>(norb_, nelea, neleb, false, true); }

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    bool operator==(const Determinants& o) const
      { return (norb_ == o.norb_ && nelea_ == o.nelea_ && neleb_ == o.neleb_ && compress_ == o.compress_); }

    // string size
    std::tuple<size_t, size_t> len_string() const { return std::make_tuple(stringa_.size(), stringb_.size()); }

    size_t lena() const { return stringa_.size(); }
    size_t lenb() const { return stringb_.size(); }

    size_t ncsfs() const;

    std::shared_ptr<Determinants> transpose() const { return std::make_shared<Determinants>(norb_, neleb_, nelea_, compress_, true); }

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
      unsigned int out = 0;
      int k = 0;
      for (int i = 0; i != norb_; ++i)
        if (bit[i]) { out += zkl(k,i, spin); ++k; }
      return out;
    }

    const std::bitset<nbit__>& stringa(int i) const { return stringa_[i]; }
    const std::bitset<nbit__>& stringb(int i) const { return stringb_[i]; }
    const std::vector<std::bitset<nbit__>>& stringa() const { return stringa_; }
    const std::vector<std::bitset<nbit__>>& stringb() const { return stringb_; }

    std::pair<std::vector<std::tuple<int, int, int>>, double> spin_adapt(const int, std::bitset<nbit__>, std::bitset<nbit__>) const;

    int nspin() const { return nelea_ - neleb_; }
    int norb()  const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }
    bool compress() const { return compress_; }

    // single index goes to normal versions (compressed based on compress_)
    const std::vector<DetMap>& phia(const int i) const { return phia_[i]; }
    const std::vector<DetMap>& phib(const int i) const { return phib_[i]; }

    // two indices goes to uncompressed versions
    const std::vector<DetMap>& phia(const int i, const int j) const { return phia_uncompressed_[i + j*norb_]; }
    const std::vector<DetMap>& phib(const int i, const int j) const { return phib_uncompressed_[i + j*norb_]; }

    const std::vector<DetMap>& phiupa(const int i) const { return phiupa_[i]; }
    const std::vector<DetMap>& phiupb(const int i) const { return phiupb_[i]; }

    const std::vector<DetMap>& phidowna(const int i) const { return phidowna_[i]; }
    const std::vector<DetMap>& phidownb(const int i) const { return phidownb_[i]; }

    std::shared_ptr<const Determinants> addalpha() const { return detaddalpha_.lock();}
    std::shared_ptr<const Determinants> remalpha() const { return detremalpha_.lock();}
    std::shared_ptr<const Determinants> addbeta() const { return detaddbeta_.lock();}
    std::shared_ptr<const Determinants> rembeta() const { return detrembeta_.lock();}

    std::shared_ptr<Determinants> addalpha() { return detaddalpha_.lock();}
    std::shared_ptr<Determinants> remalpha() { return detremalpha_.lock();}
    std::shared_ptr<Determinants> addbeta() { return detaddbeta_.lock();}
    std::shared_ptr<Determinants> rembeta() { return detrembeta_.lock();}

    template<int spin> void link(std::shared_ptr<Determinants> odet);
};


// Template function that creates the single-displacement lists (step a and b in Knowles & Handy paper).
template <int spin>
void Determinants::const_phis_(const std::vector<std::bitset<nbit__>>& string, std::vector<std::vector<DetMap>>& phi,
    std::vector<std::vector<DetMap>>& uncompressed_phi) {

  phi.clear();
  phi.resize(compress_ ? norb_*(norb_+1)/2 : norb_*norb_);
  for (auto& iter : phi ) {
    iter.reserve(string.size());
  }

  uncompressed_phi.clear();
  uncompressed_phi.resize(norb_*norb_);
  for (auto& iter : uncompressed_phi ) {
    iter.reserve(string.size());
  }

  for (auto iter = string.begin(); iter != string.end(); ++iter) {
    for (unsigned int i = 0; i != norb_; ++i) { // annihilation
      // compress_ means that we store info only for i <= j
      if ((*iter)[i] && compress_) {
        const unsigned int source = lexical<spin>(*iter);
        std::bitset<nbit__> nbit = *iter; nbit.reset(i); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          if (!(nbit[j])) {
            std::bitset<nbit__> mbit = nbit;
            mbit.set(j);
            int minij, maxij;
            std::tie(minij, maxij) = std::minmax(i,j);
            auto detmap = DetMap(lexical<spin>(mbit), sign(mbit, i, j), source);
            phi[minij+((maxij*(maxij+1))>>1)].push_back(detmap);
            uncompressed_phi[i + j*norb_].push_back(detmap);
          }
        }
      } else if ((*iter)[i]) {
        const unsigned int source = lexical<spin>(*iter);
        std::bitset<nbit__> nbit = *iter; nbit.reset(i); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          if (!(nbit[j])) {
            std::bitset<nbit__> mbit = nbit;
            mbit.set(j);
            auto detmap = DetMap(lexical<spin>(mbit), sign(mbit, i, j), source);
            phi[i + j*norb_].push_back(detmap);
            uncompressed_phi[i + j*norb_].push_back(detmap);
          }
        }
      }
    }
  }

#if 0
  // sort each vectors
  for (auto iter = phi.begin(); iter != phi.end(); ++iter) {
    std::sort(iter->begin(), iter->end());
  }
#endif
}

template<int spin> void Determinants::link(std::shared_ptr<Determinants> odet) {
  std::shared_ptr<Determinants> plusdet;
  std::shared_ptr<Determinants> det;
  if (spin==Alpha) {
    const int de = this->nelea() - odet->nelea();
    if (de == 1) std::tie(det, plusdet) = make_pair(odet, shared_from_this());
    else if (de == -1) std::tie(det, plusdet) = make_pair(shared_from_this(), odet);
    else throw std::logic_error("Determinants::link failed");
  }
  else {
    const int de = this->neleb() - odet->neleb();
    if (de == 1) std::tie(det, plusdet) = make_pair(odet, shared_from_this());
    else if (de == -1) std::tie(det, plusdet) = make_pair(shared_from_this(), odet);
    else throw std::logic_error("Determinants::link failed");
  }

  std::vector<std::vector<DetMap>> phiup;
  std::vector<std::vector<DetMap>> phidown;

  phiup.resize(norb_);
  int upsize = ( (spin==Alpha) ? plusdet->lena() : plusdet->lenb() );
  for (auto& iter : phiup) {
    iter.reserve(upsize);
  }

  phidown.resize(norb_);
  int downsize = ( (spin==Alpha) ? det->lena() : det->lenb() );
  for (auto& iter : phidown) {
    iter.reserve(downsize);
  }

  std::vector<std::bitset<nbit__>> stringplus = ( (spin==Alpha) ? plusdet->stringa() : plusdet->stringb() );
  std::vector<std::bitset<nbit__>> string = ( (spin==Alpha) ? det->stringa() : det->stringb() );

  for (auto& istring : string) {
    for (unsigned int i = 0; i != norb_; ++i) {
      if (!(istring)[i]) { // creation
        const unsigned int source = det->lexical<spin>(istring);
        std::bitset<nbit__> nbit = istring; nbit.set(i); // created.
        const unsigned int target = plusdet->lexical<spin>(nbit);
        phiup[i].push_back(DetMap(target, sign<spin>(nbit, i), source));
      }
    }
  }

  for (auto& istring : stringplus) {
    for (unsigned int i = 0; i!= norb_; ++i) {
      if (istring[i]) { // annihilation
        const unsigned int source = plusdet->lexical<spin>(istring);
        std::bitset<nbit__> nbit = istring; nbit.reset(i); //annihilated.
        const unsigned int target = det->lexical<spin>(nbit);
        phidown[i].push_back(DetMap(target, sign<spin>(nbit, i), source));
      }
    }
  }

  // finally link
  if (spin == Alpha) {
    plusdet->detremalpha_ = det;
    plusdet->phidowna_ = phidown;

    det->detaddalpha_ = plusdet;
    det->phiupa_ = phiup;
  }
  else {
    plusdet->detrembeta_ = det;
    plusdet->phidownb_ = phidown;

    det->detaddbeta_ = plusdet;
    det->phiupb_ = phiup;
  }
}

}

#endif

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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <stddef.h>
#include <memory>
#include <tuple>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <bitset>
#include <algorithm>
#include <src/util/constants.h>

namespace bagel {

struct DetMap {
  unsigned int target;
  unsigned int source;
  int sign;
  DetMap(unsigned int t, int si, unsigned int s) : target(t), source(s), sign(si) {}; 
};

// implements a determinant space
class Determinants {
  friend class Space; // TODO Is this correct?

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
    std::vector<std::bitset<nbit__> > stringa_;
    std::vector<std::bitset<nbit__> > stringb_;

    // lexical maps (Zkl)
    void const_lexical_mapping_();
    // alpha and beta string lists
    void const_string_lists_();
    // single displacement vectors Phi's
    template <int>
    void const_phis_(const std::vector<std::bitset<nbit__> >& string,
                     std::vector<std::vector<DetMap> >& target);

    // this is slow but robust implementation of bit to number converter.
    std::vector<int> bit_to_numbers(std::bitset<nbit__> bit) const {
      std::vector<int> out;
      for (int i = 0; i != norb_; ++i) if (bit[i]) out.push_back(i);
      return out;
    }

    std::bitset<nbit__> numbers_to_bit(const std::vector<int>& num) const {
      std::bitset<nbit__> out(0);
      for (auto i = num.begin(); i != num.end(); ++i) out.set(*i);
      return out;
    }

    // some utility functions
    unsigned int& zkl(int i, int j, int spin) { return zkl_[i*norb_+j+spin*nelea_*norb_]; }
    const unsigned int& zkl(int i, int j, int spin) const { return zkl_[i*norb_+j+spin*nelea_*norb_]; }

    // configuration list i^dagger j
    std::vector<std::vector<DetMap> > phia_;
    std::vector<std::vector<DetMap> > phib_;

    // configuration list i^dagger
    std::vector<std::vector<DetMap> > phiupa_;
    std::vector<std::vector<DetMap> > phiupb_;

    // configuration list i
    std::vector<std::vector<DetMap> > phidowna_;
    std::vector<std::vector<DetMap> > phidownb_;

  public:
    Determinants(const int norb, const int nelea, const int neleb, const bool compress = true);

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    // string size
    std::tuple<size_t, size_t> len_string() const { return std::make_tuple(stringa_.size(), stringb_.size()); }

    size_t lena() const { return stringa_.size(); }
    size_t lenb() const { return stringb_.size(); }

    std::string print_bit(std::bitset<nbit__> bit) const {
      std::string out;
      for (int i = 0; i != norb_; ++i) { if (bit[i]) { out += "1"; } else { out += "."; } }
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

    int sign(std::bitset<nbit__> bit, int i, int j) {
      // masking irrelevant bits
      int min, max;
      std::tie(min,max) = std::minmax(i,j);
      std::bitset<nbit__> ii(~((1 << (min+1)) - 1));
      std::bitset<nbit__> jj((1 << max) - 1); 
      bit = (bit & ii) & jj;
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

    void print(const double* const civec, const double thr) const;
    std::bitset<nbit__> stringa(int i) { return stringa_[i]; }
    std::bitset<nbit__> stringb(int i) { return stringb_[i]; }
    const std::vector<std::bitset<nbit__> >& stringa() const { return stringa_; }
    const std::vector<std::bitset<nbit__> >& stringb() const { return stringb_; }

    std::pair<std::vector<std::tuple<int, int, int> >, double> spin_adapt(const int, std::bitset<nbit__>, std::bitset<nbit__>) const;

    int nspin() const { return nelea_ - neleb_; }
    int norb()  const { return norb_; }
    int nelea() const { return nelea_; }
    int neleb() const { return neleb_; }

    const std::vector<DetMap>& phia(const int i) const { return phia_[i]; }
    const std::vector<DetMap>& phib(const int i) const { return phib_[i]; }

    const std::vector<DetMap>& phiupa(const int i) const { return phiupa_[i]; }
    const std::vector<DetMap>& phiupb(const int i) const { return phiupb_[i]; }

    const std::vector<DetMap>& phidowna(const int i) const { return phidowna_[i]; }
    const std::vector<DetMap>& phidownb(const int i) const { return phidownb_[i]; }

    std::shared_ptr<Determinants> addalpha() { return detaddalpha_.lock();}
    std::shared_ptr<Determinants> remalpha() { return detremalpha_.lock();}
    std::shared_ptr<Determinants> addbeta() { return detaddbeta_.lock();}
    std::shared_ptr<Determinants> rembeta() { return detrembeta_.lock();}
};


// Template function that creates the single-displacement lists (step a and b in Knowles & Handy paper).
template <int spin>
void Determinants::const_phis_(const std::vector<std::bitset<nbit__> >& string, std::vector<std::vector<DetMap> >& phi) {

  phi.clear();
  phi.resize(compress_ ? norb_*(norb_+1)/2 : norb_*norb_);
  for (auto iter = phi.begin(); iter != phi.end(); ++iter) {
    iter->reserve(string.size());
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
            phi[minij+((maxij*(maxij+1))>>1)].push_back(DetMap(lexical<spin>(mbit), sign(mbit, i, j), source));
          }
        }
      } else if ((*iter)[i]) {
        const unsigned int source = lexical<spin>(*iter);
        std::bitset<nbit__> nbit = *iter; nbit.reset(i); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          if (!(nbit[j])) {
            std::bitset<nbit__> mbit = nbit;
            mbit.set(j);
            phi[i+j*norb_].push_back(DetMap(lexical<spin>(mbit), sign(mbit, i, j), source));
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

}

#endif

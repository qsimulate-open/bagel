//
// Newint - Parallel electron correlation program.
// Filename: determinants.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_FCI_DETERMINANTS_H
#define __SRC_FCI_DETERMINANTS_H

#include <memory>
#include <tuple>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#ifdef USE_SSE42_INTRINSICS
#include <nmmintrin.h>
#endif

static int numofbits(unsigned int bits) {
#ifndef USE_SSE42_INTRINSICS
  bits = (bits & 0x55555555) + (bits >> 1 & 0x55555555); bits = (bits & 0x33333333) + (bits >> 2 & 0x33333333);
  bits = (bits & 0x0f0f0f0f) + (bits >> 4 & 0x0f0f0f0f); bits = (bits & 0x00ff00ff) + (bits >> 8 & 0x00ff00ff);
  return (bits & 0x0000ffff) + (bits >>16 & 0x0000ffff); // can be cheaper, but it is fine for the time being...
#else
  return _mm_popcnt_u32(bits); // i7 or later
#endif
};

// implements a determinant space
class Determinants {
  protected:
    // assuming that the number of active orbitals are the same in alpha and beta. 
    const int norb_;

    const int nelea_;
    const int neleb_;

    const bool compress_;

    // Knowles & Handy lexical mapping
    std::vector<unsigned int> zkl_; // contains zkl (Two dimenional array. See the public function).
    // string lists
    std::vector<unsigned int> stringa_;
    std::vector<unsigned int> stringb_;

    // lexical maps (Zkl)
    void const_lexical_mapping_();
    // alpha and beta string lists
    void const_string_lists_();
    // single displacement vectors Phi's
    template <int>
    void const_phis_(const std::vector<unsigned int>& string,
                     std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > >& target);

    int sign(unsigned int bit, int i, int j) {
      // masking irrelevant bits
      const unsigned int ii = ~((1 << (std::min(i,j)+1)) - 1);
      const unsigned int jj = ((1 << (std::max(i,j))) - 1); 
      bit = (bit & ii) & jj;
      return 1 - ((numofbits(bit) & 1) << 1);
    };

    // maps bit to lexical numbers.
    template <int spin> unsigned int lexical(int bit) const {
      unsigned int out = 0;
      int k = 0;
      for (int i = 0; i != norb_; ++i, bit >>= 1) 
        if (bit & 1) { out += zkl(k,i, spin); ++k; }
      return out;
    };

    // this is slow but robust implementation of bit to number converter.
    std::vector<int> bit_to_numbers(unsigned int bit) const {
      std::vector<int> out;
      for (int i = 0; bit != 0u; ++i, bit >>= 1) if (bit & 1) out.push_back(i); 
      return out;
    };

    // some utility functions
    unsigned int& zkl(int i, int j, int spin) { return zkl_[i*norb_+j+spin*nelea_*norb_]; };
    const unsigned int& zkl(int i, int j, int spin) const { return zkl_[i*norb_+j+spin*nelea_*norb_]; };

    // configuration list
    std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > > phia_;
    std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > > phib_;

  public:
    Determinants(const int norb, const int nelea, const int neleb, const bool compress = true); 
    ~Determinants() {};

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    // string size
    std::tuple<int, int> len_string() const { return std::make_tuple(stringa_.size(), stringb_.size()); }; 

    size_t lena() const { return stringa_.size(); };
    size_t lenb() const { return stringb_.size(); };

    std::string print_bit(unsigned int bit) const {
      std::string out; 
      for (int i = 0; i != norb_; ++i, bit >>=1) { if (bit&1) { out += "1"; } else { out += "."; } }
      return out;
    };
    std::string print_bit(unsigned int bit1, unsigned int bit2) const {
      std::string out; 
      for (int i = 0; i != norb_; ++i, bit1 >>=1, bit2 >>=1) {
        if (bit1&1 && bit2&1) { out += "2"; }
        else if (bit1&1) { out += "a"; }
        else if (bit2&1) { out += "b"; }
        else { out += "."; }
      }
      return out;
    };

    void print(const double* const civec, const double thr) const;
    unsigned int stringa(int i) const { return stringa_[i]; };
    unsigned int stringb(int i) const { return stringb_[i]; };
    const std::vector<unsigned int>& stringa() const { return stringa_; };
    const std::vector<unsigned int>& stringb() const { return stringb_; };

    std::pair<std::vector<std::tuple<int, int, int> >, double> spin_adapt(const int, const int, const int) const;

    int nspin() const { return nelea_ - neleb_; }; 

    const std::vector<std::tuple<unsigned int, int, unsigned int> >& phia(const int i) const { return phia_[i]; };
    const std::vector<std::tuple<unsigned int, int, unsigned int> >& phib(const int i) const { return phib_[i]; };

};


// Template function that creates the single-displacement lists (step a and b in Knowles & Handy paper).
template <int spin>
void Determinants::const_phis_(const std::vector<unsigned int>& string,
      std::vector<std::vector<std::tuple<unsigned int, int, unsigned int> > >& phi) {

  phi.clear();
  phi.resize(compress_ ? norb_*(norb_+1)/2 : norb_*norb_);
  for (auto iter = phi.begin(); iter != phi.end(); ++iter) {
    iter->reserve(string.size());
  }

  for (auto iter = string.begin(); iter != string.end(); ++iter) {
    for (unsigned int i = 0; i != norb_; ++i) { // annihilation
      const unsigned int ibit = (1 << i);
      if (ibit & *iter && compress_) {
        const unsigned int source = lexical<spin>(*iter); 
        const unsigned int nbit = (ibit^*iter); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          const unsigned int jbit = (1 << j); 
          if (!(jbit & nbit)) {
            const unsigned int mbit = jbit^nbit;
            const int minij = std::min(i,j); 
            const int maxij = std::max(i,j);
            phi[minij+((maxij*(maxij+1))>>1)].push_back(std::make_tuple(lexical<spin>(mbit), sign(mbit, i, j), source));
          }
        }
      } else if (ibit & *iter) {
        const unsigned int source = lexical<spin>(*iter); 
        const unsigned int nbit = (ibit^*iter); // annihilated.
        for (unsigned int j = 0; j != norb_; ++j) { // creation
          const unsigned int jbit = (1 << j); 
          if (!(jbit & nbit)) {
            const unsigned int mbit = jbit^nbit;
            phi[i+j*norb_].push_back(std::make_tuple(lexical<spin>(mbit), sign(mbit, i, j), source));
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
};


#endif

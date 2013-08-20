//
// BAGEL - Parallel electron correlation program.
// Filename: ras/determinants.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_RAS_DETERMINANTS_H
#define __SRC_RAS_DETERMINANTS_H

#include <memory>
#include <tuple>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <bitset>
#include <algorithm>
#include <src/util/constants.h>

namespace bagel {

class RASDeterminants {
  protected:
    struct DMap {
      const size_t source; // target is to be inferred from its location
      const size_t ij;     // displacement operators
      const int sign;      // sign of displacement

      DMap(const size_t s, const size_t ii, const int sn) : source(s), ij(ii), sign(sn) {}
    };

    std::array<int, 3> ras_;
    const int norb_;
    const int nelea_;
    const int neleb_;
    const int max_holes_;
    const int max_particles_;

    std::vector<std::vector<DMap>> phia_;
    std::vector<std::vector<DMap>> phib_;

    std::vector<std::shared_ptr<const StringSpace>> alphaspaces_;
    std::vector<std::shared_ptr<const StringSpace>> betaspaces_;

    std::vector<std::bitset<nbit__>> stringa_;
    std::vector<std::bitset<nbit__>> stringb_;

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

    const int subspace(const int nholes, const int nparticles) const { return nholes + nparticles * max_holes_; }
    const int subspace(const std::bitset<nbit__> bit) const { return nholes(bit) + nparticles(bit) * max_holes_; }

    const int nholes(const std::bitset<nbit__> bit) const { return ( norb_[0] - (bit & std::bitset<nbit__>((1ul << norb_[0]) - 1)).count() ); }
    const int nparticles(const std::bitset<nbit__> bit) const { return ( (bit & std::bitset<nbit__>(((1ul << norb_[2]) - 1) << (norb_[0] + norb_[1]))).count() ); }

    const bool allowed(const std::bitset<nbit__> bit) const { return nholes(bit) <= max_holes_ && nparticles(bit) <= max_particles_; }

  public:
    RASDeterminants(const int norb1, const int norb2, const int norb3, const int nelea, const int neleb, const int max_holes, const int max_particles, const bool mute = true);

    static const int Alpha = 0;
    static const int Beta = 1;

    bool operator==(const RASDeterminants& o) const
      { return ( nelea_ == o.nelea_ && neleb_ == o.neleb_ && max_holes_ == o.max_holes_ && max_particles_ == o.max_particles_ && ras_ == o.ras_ ); }

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

    template<int spin>
    int sign(std::bitset<nbit__> bit, int i) const {
      const std::bitset<nbit__> ii( (1 << (i)) - 1 );
      bit = bit & ii;
      return (1 - (((bit.count() + spin*nelea_) & 1 ) << 1));
    }

    int sign(std::bitset<nbit__> bit, int i, int j) const {
      // masking irrelevant bits
      int min, max;
      std::tie(min,max) = std::minmax(i,j);
      std::bitset<nbit__> ii(~((1 << (min+1)) - 1));
      std::bitset<nbit__> jj((1 << max) - 1);
      bit = (bit & ii) & jj;
      return 1 - ((bit.count() & 1) << 1);
    }

    // These access the global string lists
    const std::bitset<nbit__>& stringa(int i) const { return stringa_[i]; }
    const std::bitset<nbit__>& stringb(int i) const { return stringb_[i]; }
    const std::vector<std::bitset<nbit__>>& stringa() const { return stringa_; }
    const std::vector<std::bitset<nbit__>>& stringb() const { return stringb_; }

    const int nspin() const { return nelea_ - neleb_; }
    const int norb()  const { return norb_; }
    const int nelea() const { return nelea_; }
    const int neleb() const { return neleb_; }
    const std::array<int, 3> ras() const { return ras_; }
    const int ras(const int i) const { return ras_[i]; }

    const std::vector<DMap>& phia(const size_t ij) const { return phia_[ij]; }
    const std::vector<DMap>& phib(const size_t ij) const { return phib_[ij]; }

    template <int spin> std::shared_ptr<const StringSpace> space(const int nholes, const int nparticles)
      { return ( spin == Alpha ? alphaspaces_ : betaspaces_ )[nholes + nparticles * max_holes_]; }

    template <int spin> std::shared_ptr<const StringSpace> space(const std::bitset<nbit__> bit)
      { return space<spin>(nholes(bit), nparticles(bit)); }

    template <int spin> const size_t lexical(const std::bitset<nbit__>& bit) const { return space<spin>(bit)->lexical(bit); }

  private:
    template <int spin> void construct_phis_(const std::vector<std::bitset<nbit__>>& strings, std::vector<std::vector<DMap>>& phi);
};

template <int spin>
void construct_phis_(const std::vector<std::bitset<nbit__>>& strings, std::vector<std::vector<DMap>>& phi) {

  phi.clear();
  phi.resize( strings.size() );

  for (auto &iphi : phi) iphi.reserve(norb_ * norb_);

  auto iphi = phi.begin();
  for (auto istring = strings.begin(); istring != strings.end(); ++istring, ++iphi) {
    const std::bitset<nbit__> targetbit = *istring;
    for (int j = 0; j < norb_; ++j) {
      if ( !targetbit[j] ) continue;
      std::bitset<nbit__> intermediatebit = targetbit; intermediatebit.reset(j);
      for (int i = 0; i < norb_; ++i) {
        if ( intermediatebit[i] ) continue;
        std::bitset<nbit__> sourcebit = intermediatebit; sourcebit.set(i);
        if ( allowed(sourcebit) )
          iphi->emplace_back(lexical<spin>(sourcebit), j + i * norb_, sign(targetbit, i, j));
      }
    }
  }
}

}

#endif

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

#include <tuple>
#include <memory>
#include <unordered_map>

#include <src/ciutil/ciblock.h>
#include <src/ciutil/cistring.h>
#include <src/ciutil/cistringset.h>

namespace bagel {

using DetMapBlock = DetMapBlock_base<RASString>;
using RASBlockInfo = CIBlockInfo<RASString>;

class RASDeterminants : public std::enable_shared_from_this<RASDeterminants> {
  protected:
    std::array<int, 3> ras_;
    const int norb_;
    const int nelea_;
    const int neleb_;
    const int max_holes_;
    const int max_particles_;

    const int lenholes_; // number of combinations of holes
    const int lenparts_; // number of combinations of particles

    size_t size_;

    std::shared_ptr<const CIStringSet<RASString>> alphaspaces_;
    std::shared_ptr<const CIStringSet<RASString>> betaspaces_;

    std::vector<std::shared_ptr<const RASBlockInfo>> blockinfo_;

    std::weak_ptr<RASDeterminants> addalpha_;
    std::weak_ptr<RASDeterminants> remalpha_;
    std::weak_ptr<RASDeterminants> addbeta_;
    std::weak_ptr<RASDeterminants> rembeta_;

    std::shared_ptr<const StringMap> phia_;
    std::shared_ptr<const StringMap> phib_;

    std::vector<std::vector<DetMapBlock>> phia_ij_;
    std::vector<std::vector<DetMapBlock>> phib_ij_;

    std::vector<std::vector<DetMap>> phiupa_;
    std::vector<std::vector<DetMap>> phiupb_;

    std::vector<std::vector<DetMap>> phidowna_;
    std::vector<std::vector<DetMap>> phidownb_;

  public:
    RASDeterminants(const int norb1, const int norb2, const int norb3, const int nelea, const int neleb, const int max_holes, const int max_particles, const bool mute = false);
    RASDeterminants(std::array<int, 3> ras, const int nelea, const int neleb, const int max_holes, const int max_particles, const bool mute = false) :
      RASDeterminants(ras[0], ras[1], ras[2], nelea, neleb, max_holes, max_particles, mute) {}

    std::shared_ptr<RASDeterminants> clone(const int nelea, const int neleb) const
      { return std::make_shared<RASDeterminants>(ras_, nelea, neleb, max_holes_, max_particles_, true); }
    std::shared_ptr<RASDeterminants> transpose() const { return this->clone(neleb_, nelea_); }

    static const int Alpha = 0;
    static const int Beta = 1;

    bool operator==(const RASDeterminants& o) const
      { return ( nelea_ == o.nelea_ && neleb_ == o.neleb_ && max_holes_ == o.max_holes_ && max_particles_ == o.max_particles_ && ras_ == o.ras_ ); }


    template<int spin>
    int sign(std::bitset<nbit__> bit, int i) const {
      static_assert(nbit__ <= sizeof(unsigned long long)*8, "verify Determinant::sign (and other functions)");
      bit &= (1ull << i) - 1ull;
      return (1 - (((bit.count() + spin*nelea_) & 1 ) << 1));
    }

    int sign(std::bitset<nbit__> bit, int i, int j) const {
      // masking irrelevant bits
      int min, max;
      std::tie(min,max) = std::minmax(i,j);
      bit &= ~((1ull << (min+1)) - 1ull);
      bit &= (1ull << max) - 1ull;
      return 1 - ((bit.count() & 1) << 1);
    }

    const int nholes(const std::bitset<nbit__> bit) const { return ras_[0] - (bit & std::bitset<nbit__>((1ull << ras_[0]) - 1ull)).count(); }
    const int nparticles(const std::bitset<nbit__> bit) const { return (bit & std::bitset<nbit__>(((1ull << ras_[2]) - 1ull) << (ras_[0] + ras_[1]))).count(); }

    const bool allowed(const std::bitset<nbit__> bit) const { return nholes(bit) <= max_holes_ && nparticles(bit) <= max_particles_; }

    const bool allowed(const int nha, const int nhb, const int npa, const int npb) const
      { return ( (nha + nhb) <= max_holes_ && (npa + npb) <= max_particles_ ); }
    const bool allowed(const std::bitset<nbit__> abit, const std::bitset<nbit__> bbit) const
      { return (nholes(abit) + nholes(bbit)) <= max_holes_ && (nparticles(abit) + nparticles(bbit)) <= max_particles_; }
    const bool allowed(const std::shared_ptr<const RASString> alpha, const std::shared_ptr<const RASString> beta) const
      { return (beta->nholes() + alpha->nholes()) <= max_holes_ && (beta->nparticles() + alpha->nparticles()) <= max_particles_; }

    template <int spin>
    const std::vector<std::shared_ptr<const RASString>> allowed_spaces(std::shared_ptr<const RASString> sp) const {
      std::vector<std::shared_ptr<const RASString>> out;
      for (int jp = 0; jp + sp->nparticles() <= max_particles_; ++jp) {
        for (int ih = 0; ih + sp->nholes() <= max_holes_; ++ih) {
          std::shared_ptr<const RASString> sp = space<spin^1>(ih, jp);
          if (!sp->empty()) out.push_back(sp);
        }
      }
      return out;
    }

    // These access the global string lists
    const std::bitset<nbit__>& string_bits_a(const size_t i) const { return alphaspaces_->strings(i); }
    const std::bitset<nbit__>& string_bits_b(const size_t i) const { return betaspaces_->strings(i); }
    const std::vector<std::bitset<nbit__>>& string_bits_a() const { return alphaspaces_->strings(); }
    const std::vector<std::bitset<nbit__>>& string_bits_b() const { return betaspaces_->strings(); }

    const std::vector<std::shared_ptr<const RASBlockInfo>>& blockinfo() const { return blockinfo_; }

    const int nspin() const { return nelea_ - neleb_; }
    const int norb()  const { return norb_; }
    const int nelea() const { return nelea_; }
    const int neleb() const { return neleb_; }
    const std::array<int, 3> ras() const { return ras_; }
    const int ras(const int i) const { return ras_[i]; }

    const int max_holes() const { return max_holes_; }
    const int max_particles() const { return max_particles_; }

    const size_t lena() const { return alphaspaces_->size(); }
    const size_t lenb() const { return betaspaces_->size(); }
    const size_t size() const { return size_; }
    const int lenholes() const { return lenholes_; }
    const int lenparts() const { return lenparts_; }

    const std::vector<DetMap>& phia(const size_t target) const { return phia_->data(target); }
    const std::vector<DetMap>& phib(const size_t target) const { return phib_->data(target); }

    const std::vector<DetMapBlock>& phia_ij(const size_t ij) const { return phia_ij_[ij]; }
    const std::vector<DetMapBlock>& phib_ij(const size_t ij) const { return phib_ij_[ij]; }

    const std::vector<DetMap>& phiupa(const size_t i) const { return phiupa_[i]; }
    const std::vector<DetMap>& phiupb(const size_t i) const { return phiupb_[i]; }
    const std::vector<DetMap>& phidowna(const size_t i) const { return phidowna_[i]; }
    const std::vector<DetMap>& phidownb(const size_t i) const { return phidownb_[i]; }

    std::shared_ptr<const RASDeterminants> addalpha() const { return addalpha_.lock();}
    std::shared_ptr<const RASDeterminants> remalpha() const { return remalpha_.lock();}
    std::shared_ptr<const RASDeterminants> addbeta() const { return addbeta_.lock();}
    std::shared_ptr<const RASDeterminants> rembeta() const { return rembeta_.lock();}

    template <int spin> void link(std::shared_ptr<RASDeterminants> odet);

    template <int spin> std::shared_ptr<const RASString> space(const int nholes, const int nparticles) const {
      return (spin == Alpha ? alphaspaces_ : betaspaces_)->find_string(nholes, nparticles);
    }
    template <int spin> std::shared_ptr<const RASString> space(const std::bitset<nbit__>& bit) const {
      return (spin == Alpha ? alphaspaces_ : betaspaces_)->find_string(bit);
    }
    const std::shared_ptr<const CIStringSet<RASString>>& stringspacea() const { return alphaspaces_; }
    const std::shared_ptr<const CIStringSet<RASString>>& stringspaceb() const { return betaspaces_; }

    template <int spin> size_t lexical_zero(const std::bitset<nbit__>& bit)   const { return (spin == 0 ? alphaspaces_ : betaspaces_)->lexical_zero(bit); }
    template <int spin> size_t lexical_offset(const std::bitset<nbit__>& bit) const { return (spin == 0 ? alphaspaces_ : betaspaces_)->lexical_offset(bit); }

    std::pair<std::vector<std::tuple<std::bitset<nbit__>, std::bitset<nbit__>, int>>, double>
      spin_adapt(const int spin, const std::bitset<nbit__> alpha, const std::bitset<nbit__> beta) const;

  private:
    template <int spin> void construct_phis_(std::shared_ptr<const CIStringSet<RASString>> stringspace, std::shared_ptr<const StringMap>& phi, std::vector<std::vector<DetMapBlock>>& phi_ij);
};

template <int spin>
void RASDeterminants::construct_phis_(std::shared_ptr<const CIStringSet<RASString>> stringspace, std::shared_ptr<const StringMap>& phi, std::vector<std::vector<DetMapBlock>>& phi_ij) {

  phi = stringspace->phi();


  // old code for phi_ij TODO replace
  const int nij = (norb_ * (norb_ + 1))/2;
  const size_t stringsize = stringspace->size();
  phi_ij.clear();
  phi_ij.resize(nij);
  for (auto& iphi : phi_ij) iphi.reserve(stringsize);

  std::unordered_map<size_t, size_t> lexmap;
  for (size_t i = 0; i < stringsize; ++i)
    lexmap[(spin == 0 ? this->string_bits_a(i) : this->string_bits_b(i)).to_ullong()] = i;

  std::vector<size_t> offsets(nij, 0);

  for (auto& ispace : *stringspace) {
    size_t tindex = 0;
    for (auto istring = ispace->begin(); istring != ispace->end(); ++istring, ++tindex) {
      const std::bitset<nbit__> targetbit = *istring;
      std::vector<std::vector<DetMap>> pij;
      pij.resize( nij );
      for (int j = 0; j < norb_; ++j) {
        if ( !targetbit[j] ) continue;
        std::bitset<nbit__> intermediatebit = targetbit; intermediatebit.reset(j);
        for (int i = 0; i < norb_; ++i) {
          if ( intermediatebit[i] ) continue;
          std::bitset<nbit__> sourcebit = intermediatebit; sourcebit.set(i);
          if ( allowed(sourcebit) ) {
            const size_t source_lex = lexmap[sourcebit.to_ullong()];
            int minij, maxij;
            std::tie(minij, maxij) = std::minmax(i,j);
            pij[minij+((maxij*(maxij+1))>>1)].emplace_back(source_lex, sign(targetbit, i, j), tindex, j+i*norb_);
          }
        }
      }
      for (int i = 0; i < nij; ++i) if (pij[i].size() > 0) {
        pij[i].shrink_to_fit();
        phi_ij[i].emplace_back(offsets[i], ispace, std::move(pij[i]));
        offsets[i] += phi_ij[i].back().size();
      }
    }
  }
}

template <int spin>
void RASDeterminants::link(std::shared_ptr<RASDeterminants> odet) {
  std::shared_ptr<RASDeterminants> plusdet;
  std::shared_ptr<RASDeterminants> det;
  const int de = spin == 0 ? this->nelea() - odet->nelea() : this->neleb() - odet->neleb();
  if      (de ==  1) std::tie(det, plusdet) = std::make_pair(odet, shared_from_this());
  else if (de == -1) std::tie(det, plusdet) = std::make_pair(shared_from_this(), odet);
  else assert(false);

  std::vector<std::vector<DetMap>> phiup;
  std::vector<std::vector<DetMap>> phidown;

  phiup.resize(norb_);
  const size_t upsize = ( (spin==0) ? plusdet->lena() : plusdet->lenb() );
  for (auto& iter : phiup) {
    iter.reserve(upsize);
  }

  phidown.resize(norb_);
  const size_t downsize = ( (spin==0) ? det->lena() : det->lenb() );
  for (auto& iter : phidown) {
    iter.reserve(downsize);
  }

  std::vector<std::bitset<nbit__>> stringplus = (spin==0) ? plusdet->string_bits_a() : plusdet->string_bits_b();
  std::vector<std::bitset<nbit__>> string = (spin==0) ? det->string_bits_a() : det->string_bits_b();

  for (auto& istring : string) {
    for (unsigned int i = 0; i != norb_; ++i) {
      if (!(istring)[i]) { // creation
        const unsigned int source = det->lexical_offset<spin>(istring);
        std::bitset<nbit__> nbit = istring; nbit.set(i); // created.
        if (plusdet->allowed(nbit)) {
          const size_t target = plusdet->lexical_offset<spin>(nbit);
          phiup[i].emplace_back(target, sign<spin>(nbit, i), source, 0);
          phidown[i].emplace_back(source, sign<spin>(nbit, i), target, 0);
        }
      }
    }
  }

  // finally link
  if (spin == 0) {
    plusdet->remalpha_ = det;
    plusdet->phidowna_ = phidown;

    det->addalpha_ = plusdet;
    det->phiupa_ = phiup;
  }
  else {
    plusdet->rembeta_ = det;
    plusdet->phidownb_ = phidown;

    det->addbeta_ = plusdet;
    det->phiupb_ = phiup;
  }
}

}

#endif

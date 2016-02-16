//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/determinants.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_RAS_DETERMINANTS_H
#define __SRC_RAS_DETERMINANTS_H

#include <src/ci/ciutil/determinants_base.h>

namespace bagel {

using DetMapBlock = DetMapBlock_base<RASString>;

class RASDeterminants : public Determinants_base<RASString>,
                        public std::enable_shared_from_this<RASDeterminants> {
  protected:
    std::array<int, 3> ras_;
    const int max_holes_;
    const int max_particles_;

    std::weak_ptr<RASDeterminants> addalpha_;
    std::weak_ptr<RASDeterminants> remalpha_;
    std::weak_ptr<RASDeterminants> addbeta_;
    std::weak_ptr<RASDeterminants> rembeta_;

    std::vector<std::vector<DetMapBlock>> phia_ij_;
    std::vector<std::vector<DetMapBlock>> phib_ij_;

  public:
    RASDeterminants(const int norb1, const int norb2, const int norb3, const int nelea, const int neleb, const int max_holes, const int max_particles, const bool mute = false);
    RASDeterminants(std::array<int, 3> ras, const int nelea, const int neleb, const int max_holes, const int max_particles, const bool mute = false) :
      RASDeterminants(ras[0], ras[1], ras[2], nelea, neleb, max_holes, max_particles, mute) {}

    std::shared_ptr<RASDeterminants> clone(const int nelea, const int neleb) const
      { return std::make_shared<RASDeterminants>(ras_, nelea, neleb, max_holes_, max_particles_, true); }
    std::shared_ptr<RASDeterminants> transpose() const { return this->clone(neleb(), nelea()); }

    static const int Alpha = 0;
    static const int Beta = 1;

    int hpaddress(const int& na, const int& nb) const {
      const int N = na + nb;
      return ( (N*(N+1))/2 + nb );
    }

    int block_address(const int& nha, const int& nhb, const int& npa, const int& npb) const {
      const int lp = (max_particles()+1) * (max_particles()+2) / 2;
      return hpaddress(npa, npb) + lp * hpaddress(nha, nhb);
    }

    bool operator==(const RASDeterminants& o) const
      { return ( nelea() == o.nelea() && neleb() == o.neleb() && max_holes_ == o.max_holes_ && max_particles_ == o.max_particles_ && ras_ == o.ras_ ); }

    int nholes(const std::bitset<nbit__> bit) const { return ras_[0] - (bit & (~std::bitset<nbit__>(0ull) >> (nbit__ - ras_[0]))).count(); }
    int nparticles(const std::bitset<nbit__> bit) const { return (bit & (~(~std::bitset<nbit__>(0ull) << ras_[2]) << ras_[0] + ras_[1])).count(); }

    bool allowed(const std::bitset<nbit__> bit) const { return nholes(bit) <= max_holes_ && nparticles(bit) <= max_particles_; }

    bool allowed(const int nha, const int nhb, const int npa, const int npb) const
      { return ( (nha + nhb) <= max_holes_ && (npa + npb) <= max_particles_ ); }
    bool allowed(const std::bitset<nbit__> abit, const std::bitset<nbit__> bbit) const
      { return (nholes(abit) + nholes(bbit)) <= max_holes_ && (nparticles(abit) + nparticles(bbit)) <= max_particles_; }
    bool allowed(const std::shared_ptr<const RASString> alpha, const std::shared_ptr<const RASString> beta) const
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

    template <int spin>
    std::vector<std::shared_ptr<const CIBlockInfo<RASString>>> matching_blocks(const std::shared_ptr<const RASString>& sp) const {
      std::vector<std::shared_ptr<const CIBlockInfo<RASString>>> out;
      for (auto& binfo : blockinfo_) {
        if (!binfo->empty()) {
          if (spin==0) {
            if (sp->matches(binfo->stringsa()))
              out.push_back(binfo);
          }
          else {
            if (sp->matches(binfo->stringsb()))
              out.push_back(binfo);
          }
        }
      }
      return out;
    }

    using Determinants_base<RASString>::blockinfo;
    std::shared_ptr<const CIBlockInfo<RASString>> blockinfo(const int& nha, const int& nhb, const int& npa, const int& npb) const {
      return this->blockinfo(block_address(nha, nhb, npa, npb));
    }

    const std::array<int, 3> ras() const { return ras_; }
    int ras(const int i) const { return ras_[i]; }

    int max_holes() const { return max_holes_; }
    int max_particles() const { return max_particles_; }

    const std::vector<DetMapBlock>& phia_ij(const size_t ij) const { return phia_ij_[ij]; }
    const std::vector<DetMapBlock>& phib_ij(const size_t ij) const { return phib_ij_[ij]; }

    std::shared_ptr<const RASDeterminants> addalpha() const { return addalpha_.lock();}
    std::shared_ptr<const RASDeterminants> remalpha() const { return remalpha_.lock();}
    std::shared_ptr<const RASDeterminants> addbeta() const { return addbeta_.lock();}
    std::shared_ptr<const RASDeterminants> rembeta() const { return rembeta_.lock();}

    void set_addalpha(std::shared_ptr<RASDeterminants> o) { addalpha_ = o;}
    void set_remalpha(std::shared_ptr<RASDeterminants> o) { remalpha_ = o;}
    void set_addbeta (std::shared_ptr<RASDeterminants> o) { addbeta_ = o;}
    void set_rembeta (std::shared_ptr<RASDeterminants> o) { rembeta_ = o;}

    template <int spin> void link(std::shared_ptr<RASDeterminants> odet) { bagel::link<spin, RASString>(shared_from_this(), odet); }

    template <int spin>
    std::shared_ptr<const RASString> space(const int nholes, const int nparticles) const { return (spin == Alpha ? alphaspaces_ : betaspaces_)->find_string(nholes, nparticles); }
    template <int spin>
    std::shared_ptr<const RASString> space(const std::bitset<nbit__>& bit) const { return (spin == Alpha ? alphaspaces_ : betaspaces_)->find_string(bit); }

    std::pair<std::vector<std::tuple<std::bitset<nbit__>, std::bitset<nbit__>, int>>, double>
      spin_adapt(const int spin, const std::bitset<nbit__> alpha, const std::bitset<nbit__> beta) const;

  private:
    template <int spin> void construct_phis_(std::shared_ptr<const CIStringSet<RASString>> stringspace, std::shared_ptr<const StringMap>& phi, std::vector<std::vector<DetMapBlock>>& phi_ij);
};

template <int spin>
void RASDeterminants::construct_phis_(std::shared_ptr<const CIStringSet<RASString>> stringspace, std::shared_ptr<const StringMap>& phi, std::vector<std::vector<DetMapBlock>>& phi_ij) {

  phi = stringspace->phi();


  // old code for phi_ij TODO replace
  const int nij = (norb() * (norb() + 1))/2;
  const size_t stringsize = stringspace->size();
  phi_ij.clear();
  phi_ij.resize(nij);
  for (auto& iphi : phi_ij) iphi.reserve(stringsize);

  std::unordered_map<std::bitset<nbit__>, size_t> lexmap;
  for (size_t i = 0; i < stringsize; ++i)
    lexmap[(spin == 0 ? this->string_bits_a(i) : this->string_bits_b(i))] = i;

  std::vector<size_t> offsets(nij, 0);

  for (auto& ispace : *stringspace) {
    size_t tindex = 0;
    for (auto istring = ispace->begin(); istring != ispace->end(); ++istring, ++tindex) {
      const std::bitset<nbit__> targetbit = *istring;
      std::vector<std::vector<DetMap>> pij;
      pij.resize( nij );
      for (int j = 0; j < norb(); ++j) {
        if ( !targetbit[j] ) continue;
        std::bitset<nbit__> intermediatebit = targetbit; intermediatebit.reset(j);
        for (int i = 0; i < norb(); ++i) {
          if ( intermediatebit[i] ) continue;
          std::bitset<nbit__> sourcebit = intermediatebit; sourcebit.set(i);
          if ( allowed(sourcebit) ) {
            const size_t source_lex = lexmap[sourcebit];
            int minij, maxij;
            std::tie(minij, maxij) = std::minmax(i,j);
            pij[minij+((maxij*(maxij+1))>>1)].emplace_back(source_lex, sign(targetbit, i, j), tindex, j+i*norb());
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

}

#endif

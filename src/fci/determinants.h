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

#include <src/fci/fciblock.h>
#include <src/ciutil/cistringspace.h>

namespace bagel {

//  implements a determinant space with fixed numbers of nelea and neleb.

class Determinants : public std::enable_shared_from_this<Determinants> {
  protected:
    // vector of length 1. using unordered_map for compatibility to RAS code
    std::unordered_map<size_t, std::shared_ptr<const FCIBlockInfo>> blockinfo_;

    // global string list (TODO in principle, this is redundant)
    std::vector<std::bitset<nbit__>> string_bits_a_;
    std::vector<std::bitset<nbit__>> string_bits_b_;

    // for range check in the debug mode
    const std::shared_ptr<const FCIBlockInfo>& blockinfo(const int i) const {
      auto iter = blockinfo_.find(i);
      assert(iter != blockinfo_.end());
      return iter->second;
    }


    // TODO this should be the property of string
    bool compress_;

    /* Links to other determinant spaces accessible by one annihilation or creation operation */
    std::weak_ptr<Determinants> detaddalpha_;
    std::weak_ptr<Determinants> detaddbeta_;
    std::weak_ptr<Determinants> detremalpha_;
    std::weak_ptr<Determinants> detrembeta_;

    // single displacement vectors Phi's
    template <int>
    void const_phis_(const std::vector<std::bitset<nbit__>>& string,
                     std::vector<std::vector<DetMap>>& target, std::vector<std::vector<DetMap>>& uncompressed_target);

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

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
    }

  public:
    Determinants() { }
    Determinants(const int norb, const int nelea, const int neleb, const bool compress = true, const bool mute = false);
    Determinants(std::shared_ptr<const Determinants> o, const bool compress = true, const bool mute = false) :
      Determinants(o->norb(), o->nelea(), o->neleb(), compress, mute) {} // Shortcut to change compression of Det

    // Shortcut to make an uncompressed and muted Determinants with specified # of electrons (used for compatibility with RASDet)
    std::shared_ptr<Determinants> clone(const int nelea, const int neleb) const {
      return std::make_shared<Determinants>(norb(), nelea, neleb, false, true);
    }

    // some functions to retrieve members of FCIBlockInfo
    // block-independent members
    int norb() const { return blockinfo(0)->norb(); }
    int nelea() const { return blockinfo(0)->nelea(); }
    int neleb() const { return blockinfo(0)->neleb(); }
    int nspin() const { return nelea() - neleb(); }
    // block dependent members
    size_t lena(const size_t i = 0) const { return blockinfo(i)->lena(); }
    size_t lenb(const size_t i = 0) const { return blockinfo(i)->lenb(); }
    size_t size(const size_t i = 0) const { return blockinfo(i)->size(); }
    size_t ncsfs(const size_t i = 0) const { return blockinfo(i)->ncsfs(); }

    template<int Spin>
    int sign(const std::bitset<nbit__>& bit, const size_t pos) const {
      return blockinfo(0)->sign<Spin>(bit, pos);
    }

    static int sign(const std::bitset<nbit__>& bit, const size_t pos0, const size_t pos1) {
      return FCIBlockInfo::sign(bit, pos0, pos1);
    }

    template<int Spin>
    size_t lexical(const std::bitset<nbit__>& bit) const {
      return blockinfo(0)->lexical<Spin>(bit);
    }

    // TODO see above
    const std::bitset<nbit__>& string_bits_a(const size_t i) const { assert(i < string_bits_a_.size()); return string_bits_a_[i]; }
    const std::bitset<nbit__>& string_bits_b(const size_t i) const { assert(i < string_bits_b_.size()); return string_bits_b_[i]; }
    const std::vector<std::bitset<nbit__>>& string_bits_a() const { return string_bits_a_; }
    const std::vector<std::bitset<nbit__>>& string_bits_b() const { return string_bits_b_; }


    bool operator==(const Determinants& o) const
      { return (norb() == o.norb() && nelea() == o.nelea() && neleb() == o.neleb() && compress_ == o.compress_); }

    std::shared_ptr<Determinants> transpose() const { return std::make_shared<Determinants>(norb(), neleb(), nelea(), compress_, true); }

    std::pair<std::vector<std::tuple<int, int, int>>, double> spin_adapt(const int, std::bitset<nbit__>, std::bitset<nbit__>) const;

    bool compress() const { return compress_; }

    // single index goes to normal versions (compressed based on compress_)
    const std::vector<DetMap>& phia(const int i) const { return phia_[i]; }
    const std::vector<DetMap>& phib(const int i) const { return phib_[i]; }

    // two indices goes to uncompressed versions
    const std::vector<DetMap>& phia(const int i, const int j) const { return phia_uncompressed_[i + j*norb()]; }
    const std::vector<DetMap>& phib(const int i, const int j) const { return phib_uncompressed_[i + j*norb()]; }

    const std::vector<DetMap>& phiupa(const int i) const { return phiupa_[i]; }
    const std::vector<DetMap>& phiupb(const int i) const { return phiupb_[i]; }

    const std::vector<DetMap>& phidowna(const int i) const { return phidowna_[i]; }
    const std::vector<DetMap>& phidownb(const int i) const { return phidownb_[i]; }

    std::shared_ptr<const Determinants> addalpha() const { return detaddalpha_.lock();}
    std::shared_ptr<const Determinants> remalpha() const { return detremalpha_.lock();}
    std::shared_ptr<const Determinants> addbeta() const { return detaddbeta_.lock();}
    std::shared_ptr<const Determinants> rembeta() const { return detrembeta_.lock();}

    template<int spin> void link(std::shared_ptr<Determinants> odet);
};


// Template function that creates the single-displacement lists (step a and b in Knowles & Handy paper).
template <int spin>
void Determinants::const_phis_(const std::vector<std::bitset<nbit__>>& string, std::vector<std::vector<DetMap>>& phi,
    std::vector<std::vector<DetMap>>& uncompressed_phi) {

  phi.clear();
  phi.resize(compress_ ? norb()*(norb()+1)/2 : norb()*norb());
  for (auto& iter : phi ) {
    iter.reserve(string.size());
  }

  uncompressed_phi.clear();
  uncompressed_phi.resize(norb()*norb());
  for (auto& iter : uncompressed_phi ) {
    iter.reserve(string.size());
  }

  for (auto& istring : string) {
    for (unsigned int i = 0; i != norb(); ++i) { // annihilation
      // compress_ means that we store info only for i <= j
      if (istring[i] && compress_) {
        const unsigned int source = lexical<spin>(istring);
        std::bitset<nbit__> nbit = istring; nbit.reset(i); // annihilated.
        for (unsigned int j = 0; j != norb(); ++j) { // creation
          if (!(nbit[j])) {
            std::bitset<nbit__> mbit = nbit;
            mbit.set(j);
            int minij, maxij;
            std::tie(minij, maxij) = std::minmax(i,j);
            auto detmap = DetMap(lexical<spin>(mbit), sign(mbit, i, j), source);
            phi[minij+((maxij*(maxij+1))>>1)].push_back(detmap);
            uncompressed_phi[i + j*norb()].push_back(detmap);
          }
        }
      } else if (istring[i]) {
        const unsigned int source = lexical<spin>(istring);
        std::bitset<nbit__> nbit = istring; nbit.reset(i); // annihilated.
        for (unsigned int j = 0; j != norb(); ++j) { // creation
          if (!(nbit[j])) {
            std::bitset<nbit__> mbit = nbit;
            mbit.set(j);
            auto detmap = DetMap(lexical<spin>(mbit), sign(mbit, i, j), source);
            phi[i + j*norb()].push_back(detmap);
            uncompressed_phi[i + j*norb()].push_back(detmap);
          }
        }
      }
    }
  }

}

template<int spin> void Determinants::link(std::shared_ptr<Determinants> odet) {
  std::shared_ptr<Determinants> plusdet;
  std::shared_ptr<Determinants> det;
  if (spin==0) {
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

  phiup.resize(norb());
  int upsize = ( (spin==0) ? plusdet->lena() : plusdet->lenb() );
  for (auto& iter : phiup) {
    iter.reserve(upsize);
  }

  phidown.resize(norb());
  int downsize = ( (spin==0) ? det->lena() : det->lenb() );
  for (auto& iter : phidown) {
    iter.reserve(downsize);
  }

  std::vector<std::bitset<nbit__>> stringplus = (spin==0) ? plusdet->string_bits_a() : plusdet->string_bits_b();
  std::vector<std::bitset<nbit__>> string = (spin==0) ? det->string_bits_a() : det->string_bits_b();

  for (auto& istring : string) {
    for (unsigned int i = 0; i != norb(); ++i) {
      if (!istring[i]) { // creation
        const unsigned int source = det->lexical<spin>(istring);
        std::bitset<nbit__> nbit = istring; nbit.set(i); // created.
        const unsigned int target = plusdet->lexical<spin>(nbit);
        phiup[i].push_back(DetMap(target, sign<spin>(nbit, i), source));
        phidown[i].push_back(DetMap(source, sign<spin>(nbit, i), target));
      }
    }
  }

  // finally link
  if (spin == 0) {
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

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Determinants)

#endif

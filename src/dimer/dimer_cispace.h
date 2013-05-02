
// BAGEL - Parallel electron correlation program.
// Filename: dimer_cispace.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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



#ifndef __BAGEL_DIMER_CISPACE_H
#define __BAGEL_DIMER_CISPACE_H

#include <utility>
#include <memory>

#include <src/fci/dvec.h>
#include <src/fci/determinants.h>

namespace bagel {

class SpaceKey {
  public:
    const int S;
    const int m_s;
    const int q;

    SpaceKey(const int _s, const int _m_s, const int _q) : S(_s), m_s(_m_s), q(_q) {}

    bool operator<(const SpaceKey& o) const {
      if (S != o.S) return (S < o.S);
      else if (m_s != o.m_s) return (m_s < o.m_s);
      else return (q < o.q);
    }
    bool operator==(const SpaceKey& o) const { return ( ( (S == o.S) && (m_s == o.m_s) ) && (q == o.q) ); }

    std::string to_string() const {
      std::vector<std::string> anions = {{ "O", "A", "diA", "triA", "tetA", "pentA" }};
      std::vector<std::string> cations = {{ "O", "C", "diC", "triC", "tetC", "pentC" }};
      
      std::string out = ( q > 0 ? anions.at(q) : cations.at(-q) );
      out = out + "(2S=" + std::to_string(S) + ";2m_s=" + std::to_string(m_s) + ")";
      return out;
    }
};

class DimerCISpace {
  using SpaceMap = std::multimap<SpaceKey, std::shared_ptr<Dvec>>;
  using DMap = std::multimap<std::pair<int,int>, std::shared_ptr<Determinants>>;

  protected:
    // These are stored values of the neutral species
    std::pair<int, int> norb_;
    std::pair<int, int> nelea_;
    std::pair<int, int> neleb_;
    std::pair<int, int> nstates_;

    SpaceMap cispaceA_;
    SpaceMap cispaceB_;

    DMap detspaceA_;
    DMap detspaceB_;

  public:
    // This constructor will build the infrastructure; civecs need to be added later
    DimerCISpace(std::pair<int, int> nelea, std::pair<int, int> neleb, std::pair<int, int> norb) : norb_(norb), nelea_(nelea), neleb_(neleb) {}

    template<int unit> int norb() const { return (unit == 0 ? norb_.first : norb_.second); }
    template<int unit> int nelea() const { return (unit == 0 ? nelea_.first : nelea_.second); }
    template<int unit> int neleb() const { return (unit == 0 ? neleb_.first : neleb_.second); }
    template<int unit> int nstates() const { return (unit == 0 ? nstates_.first : nstates_.second); }

    std::pair<int, int> norb() const { return norb_; }
    std::pair<int, int> nelea() const { return nelea_; }
    std::pair<int, int> neleb() const { return neleb_; }
    std::pair<int, int> nstates() const { return nstates_; }

    template<int unit> SpaceMap& cispace() { return (unit == 0 ? cispaceA_ : cispaceB_); }

    template<int unit> std::shared_ptr<Dvec> ccvec(const int S, const int m_s, const int q) { return ccvec<unit>(SpaceKey(S,m_s,q)); }
    template<int unit> std::shared_ptr<Dvec> ccvec(SpaceKey key) {
      SpaceMap& space = (unit == 0 ? cispaceA_ : cispaceB_);
      auto iter = space.find(key);
      return (iter != space.end() ? iter->second : nullptr);
    }

    template<int unit> std::shared_ptr<Determinants> det(std::pair<const int, const int> p) { return det<unit>(p.first,p.second); }
    template<int unit> std::shared_ptr<Determinants> det(const int qa, const int qb) { 
      DMap& dets = (unit == 0 ? detspaceA_ : detspaceB_);
      auto iter = dets.find(std::make_pair(qa, qb));
      return (iter != dets.end() ? iter->second : nullptr);
    }

    template<int unit> void insert(std::shared_ptr<const Dvec> civec, const int spin = -1);
    template<int unit> std::shared_ptr<Determinants> add_det(const int qa, const int qb);
    template<int unit> std::shared_ptr<Determinants> add_det(std::pair<const int, const int> p) { return add_det<unit>(p.first,p.second); }

    void insert(std::pair<std::shared_ptr<const Dvec>, std::shared_ptr<const Dvec>> cipair) { insert<0>(cipair.first); insert<1>(cipair.second); }

    void complete();

  private:
    template<int unit> std::pair<int, int> detkey(const SpaceKey key) { return detkey<unit>(key.S, key.m_s, key.q); }
    template<int unit> std::pair<int, int> detkey(const int S, const int m_s, const int q) {
      const int nS = m_s - (nelea<unit>() - neleb<unit>());

      return std::make_pair( (nS - q)/2, -(nS + q)/2 );
    }
    template<int unit> std::pair<int, int> detkey(const int a, const int b) const {
      return std::make_pair(a - nelea<unit>(), b - neleb<unit>());
    }

    template<int unit> std::pair<int, int> detunkey(const int qa, const int qb) const {
      return std::make_pair(qa + nelea<unit>(), qb + neleb<unit>());
    }

   template<int unit> int charge(const int na, const int nb) const { return ( (nelea<unit>() + neleb<unit>()) - (na + nb) ); }
};

template<int unit> void DimerCISpace::insert(std::shared_ptr<const Dvec> civec, const int spin) {
  auto new_civec = std::make_shared<Dvec>(civec);

  const int nelea = civec->det()->nelea();
  const int neleb = civec->det()->neleb();

  const int m_s = nelea - neleb;
  const int S = (spin < 0) ? m_s : spin;
  const int Q = charge<unit>(nelea, neleb);

  // Reform Determinants object (to make sure it's the format I want)
  std::shared_ptr<Determinants> det = add_det<unit>(detkey<unit>(nelea, neleb));
  new_civec->set_det(det);

  auto& cispace = (unit == 0 ? cispaceA_ : cispaceB_);
  cispace.insert(std::make_pair(SpaceKey(S,m_s,Q), new_civec));

  const int ij = new_civec->ij();
  ((unit == 0) ? nstates_.first : nstates_.second) += ij;
}

template<int unit>
std::shared_ptr<Determinants> DimerCISpace::add_det(const int qa, const int qb) {
  const int nact = norb<unit>();

  auto& detspace = (unit == 0 ? detspaceA_ : detspaceB_);

  auto idet = detspace.find(std::make_pair(qa,qb));
  if ( idet != detspace.end()) {
    return idet->second;
  }
  else {
    int nelea, neleb;
    std::tie(nelea, neleb) = detunkey<unit>(qa,qb);
    std::shared_ptr<Determinants> det(new Determinants(nact, nelea, neleb, /*compress=*/false, /*mute=*/true));

    detspace.insert(std::make_pair(std::make_pair(qa,qb), det));

    idet = detspace.find(std::make_pair(qa+1,qb));
    if (idet != detspace.end()) det->link<0>(idet->second);

    idet = detspace.find(std::make_pair(qa-1,qb));
    if (idet != detspace.end()) det->link<0>(idet->second);

    idet = detspace.find(std::make_pair(qa,qb+1));
    if (idet != detspace.end()) det->link<1>(idet->second);

    idet = detspace.find(std::make_pair(qa,qb-1));
    if (idet != detspace.end()) det->link<1>(idet->second);

    return det;
  }
}

}

#endif

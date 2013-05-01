
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

class DimerCISpace {
  template<class T> using MMap = std::multimap<std::pair<int,int>, std::shared_ptr<T>>;

  protected:
    // These are stored values of the neutral species
    std::pair<int, int> norb_;
    std::pair<int, int> nelea_;
    std::pair<int, int> neleb_;
    std::pair<int, int> nstates_;

    MMap<Dvec> cispaceA_;
    MMap<Dvec> cispaceB_;

    MMap<Determinants> detspaceA_;
    MMap<Determinants> detspaceB_;

    bool anions_;
    bool dianions_;
    bool cations_;
    bool dications_;
    bool triplets_;
    bool quintets_;
    bool septets_;

  public:
    // This constructor will build the infrastructure; civecs need to be added later
    DimerCISpace(std::pair<int, int> nelea, std::pair<int, int> neleb, std::pair<int, int> norb) : norb_(norb), nelea_(nelea), neleb_(neleb),
      anions_(false), dianions_(false), cations_(false), dications_(false), triplets_(false), quintets_(false) {}

    template<int unit> int norb() const { return (unit == 0 ? norb_.first : norb_.second); }
    template<int unit> int nelea() const { return (unit == 0 ? nelea_.first : nelea_.second); }
    template<int unit> int neleb() const { return (unit == 0 ? neleb_.first : neleb_.second); }
    template<int unit> int nstates() const { return (unit == 0 ? nstates_.first : nstates_.second); }

    std::pair<int, int> norb() const { return norb_; }
    std::pair<int, int> nelea() const { return nelea_; }
    std::pair<int, int> neleb() const { return neleb_; }
    std::pair<int, int> nstates() const { return nstates_; }

    bool& anions() { return anions_; }
    bool& dianions() { return dianions_; }
    bool& cations() { return cations_; }
    bool& dications() { return dications_; }
    bool& triplets() { return triplets_; }
    bool& quintets() { return quintets_; }
    bool& septets() { return septets_; }

    template<int unit> std::shared_ptr<Dvec> ccvec(int qa = 0, int qb = 0) {
      MMap<Dvec>& space = (unit == 0 ? cispaceA_ : cispaceB_);
      return space.find(std::make_pair(qa,qb))->second;
    }

    template<int unit> std::shared_ptr<Determinants> det(int qa = 0, int qb = 0) { 
      MMap<Determinants>& dets = (unit == 0 ? detspaceA_ : detspaceB_);
      return dets.find(std::make_pair(qa,qb))->second; 
    }

    template<int unit> void insert(std::shared_ptr<const Dvec> civec);
    template<int unit> std::shared_ptr<Determinants> add_det(const int qa, const int qb);

    void insert(std::pair<std::shared_ptr<const Dvec>, std::shared_ptr<const Dvec>> cipair) { insert<0>(cipair.first); insert<1>(cipair.second); }

    void complete();

  private:
    template<int unit> std::pair<int, int> key(const int a, const int b) const {
      return std::make_pair(a - nelea<unit>(), b - neleb<unit>());
    };

    template<int unit> std::pair<int, int> unkey(const int qa, const int qb) const {
      return std::make_pair(qa + nelea<unit>(), qb + neleb<unit>());
    };
};

template<int unit> void DimerCISpace::insert(std::shared_ptr<const Dvec> civec) {
  auto new_civec = std::make_shared<Dvec>(civec);

  int qa, qb;
  std::tie(qa,qb) = key<unit>(civec->det()->nelea(), civec->det()->neleb());

  int Q = qa + qb;
  // Reform Determinants object (to make sure it's the format I want)
  std::shared_ptr<Determinants> det = add_det<unit>(qa,qb);
  new_civec->set_det(det);

  auto& cispace = (unit == 0 ? cispaceA_ : cispaceB_);
  cispace.insert(std::make_pair(std::make_pair(qa,qb), new_civec));

  const int ij = new_civec->ij();
  if (unit == 0) nstates_ = std::make_pair(nstates_.first + ij, nstates_.second);
  else nstates_ = std::make_pair(nstates_.first, nstates_.second + ij);
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
    std::tie(nelea, neleb) = unkey<unit>(qa,qb);
    auto det = std::make_shared<Determinants>(nact, nelea, neleb, /*compress=*/false, /*mute=*/true);

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

//
// BAGEL - Parallel electron correlation program.
// Filename: space.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu >
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

#include <stdexcept>
#include <src/zfci/relspace.h>
#include <src/math/comb.h>
#include <src/util/combination.hpp>

using namespace std;
using namespace bagel;

RelSpace::RelSpace(const int norb, const int nelea, const int neleb, const bool mute, const bool linkup)
  : Space_base(norb, nelea, neleb, mute), linkup_(linkup) {

  common_init();
}


void RelSpace::common_init() {
  assert(nelea_ == neleb_);

  if (!mute_) cout << " Constructing " << nelea_+neleb_+1 << " determinant spaces" << endl;
  const int nele = nelea_+neleb_;
  for (int i = 0; i != norb_; ++i) {
    if (nele-i >= 0 && nele-i <= norb_) {
      if (!mute_) cout << " Constructing space of all determinants with " << nelea_+i << " " << neleb_+i << " "  << endl << endl;
      { 
        auto tmpdet = make_shared<Determinants>(norb_, i, nele-i, false, mute_);
        detmap_.insert(pair<int,shared_ptr<Determinants>>(key_(i, nele-i), tmpdet));
      }
      if (linkup_) {
        if (i+1 <= norb_) {
          auto tmpdet = make_shared<Determinants>(norb_, i+1, nele-i, false, mute_);
          detmap_.insert(pair<int,shared_ptr<Determinants>>(key_(i+1, nele-i), tmpdet));
        }
        if (nele-i+1 <= norb_) {
          auto tmpdet2 = make_shared<Determinants>(norb_, i, nele-i+1, false, mute_);
          detmap_.insert(pair<int,shared_ptr<Determinants>>(key_(i, nele-i+1), tmpdet2));
        }
      }
    }
  }

  if (!mute_) cout << " Space is made up of " << detmap_.size() << " determinants." << endl;

  if (linkup_) {
    if (!mute_) cout << "  o forming alpha links" << endl;

    int nlinks = 0;
    for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
      int na = idet->second->nelea(); int nb = idet->second->neleb();
      auto jdet = detmap_.find(key_(na+1,nb));
      if (jdet != detmap_.end()) {
        idet->second->link<0>(jdet->second);
        ++nlinks;
      }
    }
    if (!mute_) cout << "    - " << nlinks << " links formed" << endl;

    if (!mute_) cout << "  o forming beta links" << endl;

    nlinks = 0;
    for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
      int na = idet->second->nelea(); int nb = idet->second->neleb();
      auto jdet = detmap_.find(key_(na,nb+1));
      if (jdet != detmap_.end()) {
        idet->second->link<1>(jdet->second);
        ++nlinks;
      }
    }
    if (!mute_) cout << "    - " << nlinks << " links formed" << endl;
  }
}

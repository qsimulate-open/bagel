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

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::RelSpace)

using namespace std;
using namespace bagel;

RelSpace::RelSpace(const int norb, const int nele, const bool mute, const bool linkup) {

  for (int i = 0; i <= norb; ++i) {
    if (nele-i >= 0 && nele-i <= norb) {
      if (!mute) cout << " Constructing space of all determinants with " << i << " " << nele-i << " "  << endl << endl;
      {
        auto tmpdet = make_shared<Determinants>(norb, i, nele-i, false, mute);
        detmap_.insert(pair<int,shared_ptr<Determinants>>(key_(i, nele-i), tmpdet));
      }
      if (linkup) {
        if (i+1 <= norb) {
          auto tmpdet = make_shared<Determinants>(norb, i+1, nele-i, false, mute);
          detmap_.insert(pair<int,shared_ptr<Determinants>>(key_(i+1, nele-i), tmpdet));
        }
        if (nele-i+1 <= norb) {
          auto tmpdet2 = make_shared<Determinants>(norb, i, nele-i+1, false, mute);
          detmap_.insert(pair<int,shared_ptr<Determinants>>(key_(i, nele-i+1), tmpdet2));
        }
      }
    }
  }

  if (!mute) cout << " Space is made up of " << detmap_.size() << " determinants." << endl;

  if (linkup) {
    if (!mute) cout << "  o forming alpha links" << endl;

    int nlinks = 0;
    for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
      int na = idet->second->nelea(); int nb = idet->second->neleb();
      auto jdet = detmap_.find(key_(na+1,nb));
      if (jdet != detmap_.end()) {
        idet->second->link<0>(jdet->second);
        ++nlinks;
      }
    }
    if (!mute) cout << "    - " << nlinks << " links formed" << endl;

    if (!mute) cout << "  o forming beta links" << endl;

    nlinks = 0;
    for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
      int na = idet->second->nelea(); int nb = idet->second->neleb();
      auto jdet = detmap_.find(key_(na,nb+1));
      if (jdet != detmap_.end()) {
        idet->second->link<1>(jdet->second);
        ++nlinks;
      }
    }
    if (!mute) cout << "    - " << nlinks << " links formed" << endl;
  }
}

//
// BAGEL - Parallel electron correlation program.
// Filename: space.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker < shane.parker@u.northwestern.edu >
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
#include <src/fci/determinants.h>
#include <src/fci/space.h>
#include <src/math/comb.h>
#include <src/util/combination.hpp>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::HZSpace)

using namespace std;
using namespace bagel;

HZSpace::HZSpace(const int norb, const int nelea, const int neleb, const bool compress, const bool mute) {

  if (!mute) cout << " Constructing space of all determinants that can formed by removing 1 electron from " << nelea
                  << " alpha and " << neleb << " beta electrons." << endl << endl;

  assert(neleb >= 1 && nelea >= 1);
  for(int i = -1; i <= 0; ++i ) {
    for(int j = -1; j <= 0; ++j) {
      detmap_.insert(make_pair(key_(nelea+i,neleb+j), make_shared<Determinants>(norb, nelea + i, neleb + j, compress, /*mute_=*/true)));
    }
  }
  if (!mute) {
     cout << " Space is made up of " << detmap_.size() << " determinants." << endl;
     cout << "  o forming alpha links" << endl;
  }

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

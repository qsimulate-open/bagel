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

using namespace std;
using namespace bagel;

Space::Space(const int _norb, const int _nelea, const int _neleb, const int _M, const bool _compress, const bool _mute)
  : Space_Base(_norb, _nelea, _neleb, _mute), M_(_M), compress_(_compress) {

  common_init();
}

Space::Space(shared_ptr<const Determinants> det, int _M, const bool _compress, const bool _mute) :
  Space_Base(det, _mute), M_(_M), compress_(_compress) {

  common_init();
}

void Space::common_init() {
  if (!mute_) cout << " Constructing space of all determinants that can formed by removing "
                  << M_ << " electrons from " << nelea_
                  << " alpha and " << neleb_ << " beta electrons." << endl << endl;
  for(int i = -M_; i <= 0; ++i ) {
    for(int j = -M_; j <= 0; ++j) {
      auto tmpdet = make_shared<Determinants>(norb_, nelea_ + i, neleb_ + j, compress_, /*mute_=*/true);
      detmap_.insert(pair<int,shared_ptr<Determinants>>(key_(i,j), tmpdet));
    }
  }
  if (!mute_) cout << " Space is made up of " << detmap_.size() << " determinants." << endl;
  if (!mute_) cout << "  o forming alpha links" << endl;

  int nlinks = 0;
  for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
    int na = idet->second->nelea(); int nb = idet->second->neleb();
    auto jdet = detmap_.find(key_(na-nelea_+1,nb-neleb_));
    if(jdet==detmap_.end()) continue;
    else {
      idet->second->link<0>(jdet->second);
      ++nlinks;
    }
  }
  if (!mute_) cout << "    - " << nlinks << " links formed" << endl;

  if (!mute_) cout << "  o forming beta links" << endl;

  nlinks = 0;
  for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
    int na = idet->second->nelea(); int nb = idet->second->neleb();
    auto jdet = detmap_.find(key_(na-nelea_,nb-neleb_+1));
    if(jdet==detmap_.end()) continue;
    else {
      idet->second->link<1>(jdet->second);
      ++nlinks;
    }
  }
  if (!mute_) cout << "    - " << nlinks << " links formed" << endl;
}

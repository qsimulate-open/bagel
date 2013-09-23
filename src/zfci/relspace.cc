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

RelSpace::RelSpace(const int norb, const int nelea, const int neleb, const bool mute)
  : Space_Base(norb, nelea, neleb, mute) {

  common_init();
}

RelSpace::RelSpace(shared_ptr<const Determinants> _det, const bool _mute) :
  Space_Base(_det, _mute) {

  common_init();
}

void RelSpace::common_init() {
  if (!mute_) cout << " **Constructing " << 2*(norb_ - nelec()/2)+1 << " determinant spaces by Kramers index" << endl;
  for (int M_ = -(norb_-nelec()/2) ; M_ <= norb_ - nelec()/2; ++M_) {
    if (!mute_) cout << " **Constructing space of all determinants with Kramers index "
                     << kramers(M_) << endl << endl;
    auto tmpdet = make_shared<Determinants>(norb_, nelea_ + M_, neleb_ - M_, false, mute_);
    detmap_.insert(pair<int,shared_ptr<Determinants>>(kramers(M_), tmpdet));
  }
    if (!mute_) cout << " **Total space is made up of " << detmap_.size() << " determinants." << endl;
#if 0
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
#endif
}

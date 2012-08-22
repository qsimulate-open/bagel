//
// Newint - Parallel electron correlation program.
// Filename: space.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker < shane.parker@u.northwestern.edu >
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <stdexcept>
#include <src/newfci/determinants.h>
#include <src/newfci/space.h>
#include <src/util/comb.h>
#include <src/util/combination.hpp>

using namespace std;

Space::Space(const int _norb, const int _nelea, const int _neleb, const int _M, const bool _compress)
  : norb_(_norb), nelea_(_nelea), neleb_(_neleb), M_(_M), compress_(_compress) {

  const bool mute = !compress_;

  if (!mute) cout << " Constructing space of all determinants that can formed by adding or removing " 
                  << M_ << " electrons from " << nelea_ 
                  << " alpha and " << neleb_ << " beta electrons." << endl << endl;
  for(int i = -M_; i <= M_; ++i ) {
    for(int j = -M_; j <= M_; ++j) {
      if ( ::abs(i+j) > M_ ) continue;
      else {
        shared_ptr<NewDeterminants> tmpdet(new NewDeterminants(norb_, nelea_ + i, neleb_ + j, compress_));
        detmap_.insert(pair<int,shared_ptr<NewDeterminants> >(key_(nelea_ + i,neleb_ + j), tmpdet));
      }
    }
  }
  if (!mute) cout << " Space is made up of " << detmap_.size() << " determinants." << endl;
  if (!mute) cout << "  o forming alpha links" << endl;

  int nlinks = 0;
  for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
    int na = idet->second->nelea(); int nb = idet->second->neleb();
    auto jdet = detmap_.find(key_(na+1,nb));
    if(jdet==detmap_.end()) continue;
    else {
      form_link_<0>(idet->second, jdet->second);
      ++nlinks;
    }
  }
  if (!mute) cout << "    - " << nlinks << " links formed" << endl;

  if (!mute) cout << "  o forming beta links" << endl;

  nlinks = 0;
  for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
    int na = idet->second->nelea(); int nb = idet->second->neleb();
    auto jdet = detmap_.find(key_(na,nb+1));
    if(jdet==detmap_.end()) continue;
    else {
      form_link_<1>(idet->second, jdet->second);
      ++nlinks;
    }
  }
  if (!mute) cout << "    - " << nlinks << " links formed" << endl;
}


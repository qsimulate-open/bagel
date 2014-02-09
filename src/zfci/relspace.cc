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

  map<int, shared_ptr<const FCIString>> s;
  list<shared_ptr<const FCIString>> lista;
  for (int i = 0; i <= norb; ++i)
    if ((nele-i >= 0 && nele-i <= norb) || (linkup && nele-i+1 >= 0 && nele-i+1 <= norb)) {
      s[i] = make_shared<FCIString>(i, norb);
      lista.push_back(s[i]);
    }

  auto space = make_shared<CIStringSpace<FCIString>>(lista);
  space->build_linkage();

  // make Nele determinants
  for (int i = 0; i <= norb; ++i) {
    if (nele-i >= 0 && nele-i <= norb) {
      detmap_.insert(make_pair(key_(i, nele-i), make_shared<Determinants>(s.at(i), s.at(nele-i), false, mute)));

      if (linkup) {
        if (i+1 <= norb)
          detmap_.insert(make_pair(key_(i+1, nele-i), make_shared<Determinants>(s.at(i+1), s.at(nele-i), false, mute)));
        if (nele-i+1 <= norb)
          detmap_.insert(make_pair(key_(i, nele-i+1), make_shared<Determinants>(s.at(i), s.at(nele-i+1), false, mute)));
      }
    }
  }

  if (!mute) cout << " Space is made up of " << detmap_.size() << " determinants." << endl;

  if (linkup) {
    if (!mute) cout << "  o forming links" << endl;

    for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
      auto jdet = idet; ++jdet;
      for ( ; jdet != detmap_.end(); ++jdet)
        idet->second->link(jdet->second, space, space);
    }
  }
}

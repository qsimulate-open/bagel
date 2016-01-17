//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relspace.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu >
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <stdexcept>
#include <src/ci/zfci/relspace.h>
#include <src/util/math/comb.h>
#include <src/util/combination.hpp>

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::RelSpace)

using namespace std;
using namespace bagel;

RelSpace::RelSpace(const int norb, const int nele, const bool mute, const bool linkup) : nele_(nele) {

  using FCIStringSet = CIStringSet<FCIString>;

  map<int, shared_ptr<const FCIStringSet>> s;
  list<shared_ptr<const FCIStringSet>> lista;
  for (int i = 0; i <= norb; ++i)
    if ((nele-i >= 0 && nele-i <= norb) || (linkup && nele-i+1 >= 0 && nele-i+1 <= norb)) {
      auto tmp = make_shared<FCIString>(i, norb);
      s[i] = make_shared<FCIStringSet>(list<shared_ptr<const FCIString>>{tmp});
      lista.push_back(s[i]);
    }

  auto space = make_shared<CIStringSpace<FCIStringSet>>(lista);
  space->build_linkage();

  // make Nele determinants
  for (int i = 0; i <= norb; ++i) {
    if (nele-i >= 0 && nele-i <= norb) {
      detmap_.emplace(make_pair(i, nele-i), make_shared<Determinants>(s.at(i), s.at(nele-i), false, mute));

      if (linkup) {
        if (i+1 <= norb)
          detmap_.emplace(make_pair(i+1, nele-i), make_shared<Determinants>(s.at(i+1), s.at(nele-i), false, mute));
        if (nele-i+1 <= norb)
          detmap_.emplace(make_pair(i, nele-i+1), make_shared<Determinants>(s.at(i), s.at(nele-i+1), false, mute));
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

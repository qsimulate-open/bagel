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

  list<shared_ptr<const FCIString>> lista{make_shared<FCIString>(nelea-1, norb), make_shared<FCIString>(nelea, norb)};
  auto spacea = make_shared<CIStringSpace<FCIString>>(lista);
  spacea->build_linkage();

  list<shared_ptr<const FCIString>> listb{make_shared<FCIString>(neleb-1, norb), make_shared<FCIString>(neleb, norb)};
  auto spaceb = make_shared<CIStringSpace<FCIString>>(listb);
  spaceb->build_linkage();

  if (!mute) {
    cout << " Space is made up of " << detmap_.size() << " determinants." << endl;
    cout << "  o forming links" << endl;
  }

  // build determinants
  for (auto& a : lista)
    for (auto& b : listb)
      detmap_.insert(make_pair(key_(a->nele(), b->nele()), make_shared<Determinants>(a, b, compress, true)));

  for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
    auto jdet = idet; ++jdet;
    for( ; jdet != detmap_.end(); ++jdet)
      idet->second->link(jdet->second, spacea, spaceb);
  }

}

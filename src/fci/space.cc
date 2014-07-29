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
  auto na0 = make_shared<FCIString>(nelea  , norb);
  auto na1 = make_shared<FCIString>(nelea-1, norb);
  auto nb0 = make_shared<FCIString>(neleb  , norb);
  auto nb1 = make_shared<FCIString>(neleb-1, norb);

  using FCIStringSet = CIStringSet<FCIString>;

  list<shared_ptr<const FCIStringSet>> lista = { make_shared<FCIStringSet>(list<shared_ptr<const FCIString>>{na0}),
                                                 make_shared<FCIStringSet>(list<shared_ptr<const FCIString>>{na1}) };
  list<shared_ptr<const FCIStringSet>> listb = { make_shared<FCIStringSet>(list<shared_ptr<const FCIString>>{nb0}),
                                                 make_shared<FCIStringSet>(list<shared_ptr<const FCIString>>{nb1}) };

  spacea_ = make_shared<CIStringSpace<FCIStringSet>>(lista);
  spacea_->build_linkage();

  spaceb_ = make_shared<CIStringSpace<FCIStringSet>>(listb);
  spaceb_->build_linkage();

  if (!mute) {
    cout << " Space is made up of " << detmap_.size() << " determinants." << endl;
    cout << "  o forming links" << endl;
  }

  // build determinants
  for (auto& a : lista)
    for (auto& b : listb)
      detmap_.emplace(make_pair(a->nele(), b->nele()), make_shared<Determinants>(a, b, compress, true));

  // link determinants
  link();
}


void HZSpace::link() {
  for(auto idet = detmap_.begin(); idet != detmap_.end(); ++idet) {
    auto jdet = idet; ++jdet;
    for( ; jdet != detmap_.end(); ++jdet)
      idet->second->link(jdet->second, spacea_, spaceb_);
  }
}

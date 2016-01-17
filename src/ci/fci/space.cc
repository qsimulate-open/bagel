//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: space.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker < shane.parker@u.northwestern.edu >
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
#include <src/ci/fci/determinants.h>
#include <src/ci/fci/space.h>
#include <src/util/math/comb.h>
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

  // build determinants
  for (auto& a : lista)
    for (auto& b : listb)
      detmap_.emplace(make_pair(a->nele(), b->nele()), make_shared<Determinants>(a, b, compress, true));

  if (!mute) {
    cout << " Space is made up of " << detmap_.size() << " determinants." << endl;
    cout << "  o forming links" << endl;
  }

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

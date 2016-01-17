//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: phi_k_lists.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/dmrg/phi_k_lists.h>

using namespace std;
using namespace bagel;

PhiKLists::PhiKLists(shared_ptr<const CIStringSet<RASString>> source_stringspace, shared_ptr<const CIStringSet<RASString>> target_stringspace) {
  assert(source_stringspace->norb() == target_stringspace->norb());
  const int norb = source_stringspace->norb();
  const int nele_target = target_stringspace->nele();

  for (int k = 0; k < norb; ++k) {
    unordered_map<int, vector<PhiK>> phimap;

    for (auto& source_space : *source_stringspace) {
      vector<PhiK> iphilist;
      for (size_t ia = 0; ia < source_space->size(); ++ia) {
        const bitset<nbit__> sbit = source_stringspace->strings(ia + source_space->offset());
        const bitset<nbit__> tbit = sbit ^ bitset<nbit__>(1 << k);

        // counting nelea dictates whether the target bit belongs to the right set of determinants
        if (tbit.count() == nele_target && target_stringspace->allowed(tbit)) {
          // assume these are alpha bits for the purposes of phase
          const int phase = sign(tbit, k);
          const size_t target_lex = target_stringspace->lexical_zero(tbit);
          iphilist.emplace_back(ia, target_lex, phase, target_stringspace->find_string(tbit).get());
        }
      }
      phimap.emplace(source_space->tag(), move(iphilist));
    }
    data_.emplace(k, move(phimap));
  }
}

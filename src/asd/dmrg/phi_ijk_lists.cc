//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: phi_ijk_lists.cc
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

#include <src/asd/dmrg/phi_ijk_lists.h>

using namespace std;
using namespace bagel;

PhiIJKLists::PhiIJKLists(shared_ptr<const CIStringSet<RASString>> source_stringspace,
            shared_ptr<const CIStringSet<RASString>> target_stringspace, const bool conjugate) {
  data_.reserve(target_stringspace->size());
  const size_t tlen = target_stringspace->size();
  const int norb = source_stringspace->norb();
  assert(norb == target_stringspace->norb());

  if (conjugate) {
    for (size_t it = 0; it < tlen; ++it) {
      vector<PhiIJK> phi;
      const bitset<nbit__> tbit = target_stringspace->strings(it);
      for (int k = 0; k < norb; ++k) {
        if (!tbit[k]) continue;
        const int kphase = sign(tbit, k);
        for (int j = 0; j < k; ++j) {
          if (!tbit[j]) continue;
          const bitset<nbit__> tmpbit = tbit ^ ((bitset<nbit__>(1) << k) | (bitset<nbit__>(1) << j));
          for (int i = 0; i < norb; ++i) {
            if (tmpbit[i]) continue;
            const bitset<nbit__> sbit = tmpbit ^ (bitset<nbit__>(1) << i);
            if (source_stringspace->allowed(sbit)) {
              const int ijphase = sign(sbit, i, j);
              const size_t source_lex = source_stringspace->lexical_offset(sbit);
              phi.emplace_back(source_lex, ijphase*kphase, i, j, k);
            }
          }
        }
      }
      data_.emplace_back(move(phi));
    }
  } else {
    for (size_t it = 0; it < tlen; ++it) {
      vector<PhiIJK> phi;
      const bitset<nbit__> tbit = target_stringspace->strings(it);
      for (int i = 0; i < norb; ++i) {
        if (!tbit[i]) continue;
        const bitset<nbit__> tmpbit = tbit ^ (bitset<nbit__>(1) << i);
        for (int j = 0; j < norb; ++j) {
          if (tmpbit[j]) continue;
          const bitset<nbit__>tmp2bit = tmpbit ^ (bitset<nbit__>(1) << j);
          const int ijphase = sign(tbit, i, j);
          for (int k = j+1; k < norb; ++k) {
            if (tmpbit[k]) continue;
            const bitset<nbit__> sbit = tmp2bit ^ (bitset<nbit__>(1) << k);
            if (source_stringspace->allowed(sbit)) {
              const int kphase = sign(sbit, k);
              const size_t source_lex = source_stringspace->lexical_offset(sbit);
              phi.emplace_back(source_lex, ijphase*kphase, i, j, k);
            }
          }
        }
      }
      data_.emplace_back(move(phi));
    }
  }
}

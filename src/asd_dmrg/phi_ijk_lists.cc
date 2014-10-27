//
// BAGEL - Parallel electron correlation program.
// Filename: phi_ijk_lists.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd_dmrg/phi_ijk_lists.h>

using namespace std;
using namespace bagel;

PhiIJKLists::PhiIJKLists(shared_ptr<const CIStringSet<RASString>> source_stringspace,
            shared_ptr<const CIStringSet<RASString>> target_stringspace, const bool conjugate)
{
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
  }
  else {
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

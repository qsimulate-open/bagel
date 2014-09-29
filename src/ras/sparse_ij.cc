//
// BAGEL - Parallel electron correlation program.
// Filename: ras/sparse_ij.cc
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

#include <src/ras/sparse_ij.h>

using namespace std;
using namespace bagel;

Sparse_IJ::Sparse_IJ(shared_ptr<const CIStringSet<RASString>> source_stringspace, shared_ptr<const CIStringSet<RASString>> target_stringspace) {
  assert(source_stringspace->nele()==target_stringspace->nele());
  assert(source_stringspace->norb()==target_stringspace->norb());

  const int norb = source_stringspace->norb();

  for (auto source_space : *source_stringspace) {
    for (auto target_space : *target_stringspace) {
      map<pair<int, int>, double> sparse_coords;
      map<pair<int, int>, vector<SparseIJKey>> sparse_keys;

      for (size_t itar = 0; itar < target_space->size(); ++itar) {
        const bitset<nbit__> tbit = target_space->strings(itar);
        for (int i = 0; i < norb; ++i) {
          if (!tbit[i]) continue;
          const bitset<nbit__> tmpbit = tbit ^ bitset<nbit__>(1 << i);
          for (int j = 0; j < norb; ++j) {
            if (tmpbit[j]) continue;
            const bitset<nbit__> sbit = tmpbit ^ bitset<nbit__>(1 << j);
            if (source_space->contains(sbit)) {
              const size_t isrc = source_space->lexical_zero(sbit);
              sparse_coords[{itar,isrc}] = 0.0;
              sparse_keys[{itar,isrc}].emplace_back(i, j, sign(sbit, i, j), nullptr);
            }
          }
        }
      }

      if (sparse_coords.size() > 0) {
        auto sparse = make_shared<SparseMatrix>(target_space->size(), source_space->size(), sparse_coords);
        double* data = sparse->data();
        vector<SparseIJKey> sp;
        sp.reserve(accumulate(sparse_keys.begin(), sparse_keys.end(), 0, [] (int x, const pair<pair<int,int>, vector<SparseIJKey>>& i) { return x + i.second.size(); }));

        for (auto i : sparse_keys) {
          for (auto d : i.second)
            sp.emplace_back(d.i, d.j, d.sign, data);
          ++data;
        }
        data_.emplace(make_pair(target_space->tag(), source_space->tag()), make_tuple(sparse, move(sp)));
      }
      else {
        data_.emplace(make_pair(target_space->tag(), source_space->tag()), make_tuple(nullptr, vector<SparseIJKey>()));
      }
    }
  }
}

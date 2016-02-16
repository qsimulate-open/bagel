//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/sparse_ij.cc
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

#include <set>
#include <src/ci/ras/sparse_ij.h>

using namespace std;
using namespace bagel;

Sparse_IJ::Sparse_IJ(shared_ptr<const CIStringSet<RASString>> source_stringspace, shared_ptr<const CIStringSet<RASString>> target_stringspace) {
  assert(source_stringspace->nele()==target_stringspace->nele());
  assert(source_stringspace->norb()==target_stringspace->norb());

  const int norb = source_stringspace->norb();
  unordered_map<bitset<nbit__>, size_t> lexmap;
  for (auto source_space : *source_stringspace)
    for (size_t i = 0; i < source_space->size(); ++i)
      lexmap[source_space->strings(i)] = i;

  const size_t slen = source_stringspace->size();

  for (auto& target_space : *target_stringspace) {
    const int nsource = source_stringspace->nspaces();
    vector<vector<tuple<size_t, SparseIJKey>>> sparse_keys_vec(nsource);
    vector<set<size_t>> sparse_set_vec(nsource);

    for (size_t itar = 0; itar < target_space->size(); ++itar) {
      const bitset<nbit__> tbit = target_space->strings(itar);
      for (int i = 0; i < norb; ++i) {
        if (!tbit[i]) continue;
        const bitset<nbit__> tmpbit = tbit ^ bitset<nbit__>(1 << i);
        for (int j = 0; j < norb; ++j) {
          if (tmpbit[j]) continue;
          const bitset<nbit__> sbit = tmpbit ^ bitset<nbit__>(1 << j);
          int isource_space = 0;
          for (auto& source_space : *source_stringspace) {
            if (source_space->contains(sbit)) {
              const size_t isrc = lexmap[sbit];
              sparse_set_vec[isource_space].emplace(isrc + slen * itar);
              sparse_keys_vec[isource_space].emplace_back(isrc + slen * itar, SparseIJKey(i, j, sign(sbit, i, j), nullptr));
              break;
            }
            ++isource_space;
          }
        }
      }
    }

    int isource_space = 0;
    for (auto& source_space : *source_stringspace) {
      vector<tuple<size_t, SparseIJKey>>& sparse_keys = sparse_keys_vec[isource_space];
      if (sparse_keys.size() > 0) {
        const set<size_t>& sparse_set = sparse_set_vec[isource_space];
        vector<tuple<int, int, double>> sparse_coords;
        sparse_coords.reserve(sparse_set.size());
        for (auto& i : sparse_set)
          sparse_coords.emplace_back(i/slen, i%slen, 0.0);
        auto sparse = make_shared<SparseMatrix>(target_space->size(), source_space->size(), sparse_coords);
        double* data = sparse->data();

        sort(sparse_keys.begin(), sparse_keys.end(), [] (const tuple<size_t, SparseIJKey>& a, const tuple<size_t, SparseIJKey>& b) {
          return get<0>(a) < get<0>(b);
        });

        vector<SparseIJKey> sp;
        sp.reserve(sparse_keys.size());
        size_t last = get<0>(sparse_keys.front());
        for (auto dd : sparse_keys) {
          const SparseIJKey& d = get<1>(dd);
          if (last != get<0>(dd))
            ++data;
          sp.emplace_back(d.i, d.j, d.sign, data);
          last = get<0>(dd);
        }
        data_.emplace(make_pair(target_space->tag(), source_space->tag()), make_tuple(sparse, move(sp)));
      }
      else {
        data_.emplace(make_pair(target_space->tag(), source_space->tag()), make_tuple(nullptr, vector<SparseIJKey>()));
      }
      ++isource_space;
    }
  }
}

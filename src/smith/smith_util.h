//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: smith_util.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_SMITH_SMITH_UTIL_H
#define __SRC_SMITH_SMITH_UTIL_H

#include <src/smith/tensor.h>
#include <src/util/kramers.h>

namespace bagel {
namespace SMITH {

template<int N, typename DataType>
static void fill_block(std::shared_ptr<TATensor<DataType,N>> target, std::shared_ptr<const btas::TensorN<DataType,N>> input,
                       const std::vector<int>& inpoffsets_rev) {
  assert(input->range().ordinal().contiguous());
  assert(N == input->range().rank());
  const int rank = N;

  auto prod = [](const size_t n, const Index& i) { return n*i.size(); };

  // first index corresponds to the fastest in the data
  const std::vector<IndexRange> ranges_rev = target->indexrange();
  // last index corresponds to the fastest
  const std::vector<IndexRange> ranges(ranges_rev.rbegin(), ranges_rev.rend());
  // offset has to be reversed as well. last index corresponds to the fastest
  const std::vector<int> inpoffsets(inpoffsets_rev.rbegin(), inpoffsets_rev.rend());

  // The last element corresponds to the fastest index
  std::vector<std::map<size_t,Index>> keymap;
  for (auto it = ranges.begin(); it != ranges.end(); ++it) {
    std::map<size_t,Index> key;
    size_t off = 0lu;
    for (auto& j : *it) {
      key.emplace(off, j);
      off += j.size();
    }
    keymap.push_back(key);
  }

  // loop over tiles of TATensor
  for (auto it = target->begin(); it != target->end(); ++it) {
    const TiledArray::Range range = target->trange().make_tile_range(it.ordinal());
    auto lo = range.lobound();
    assert(lo.size() == N);
    // this vector<Index> is what we use in SMITH tensor
    // the last Index correponds to the fastest index
    std::vector<Index> indices(lo.size());
    {
      auto iter = indices.begin();
      auto key = keymap.begin();
      for (auto jt = lo.begin(); jt != lo.end(); ++jt, ++key, ++iter)
        *iter = key->at(*jt);
    }

    std::vector<size_t> stride;
    for (auto i = indices.begin(); i != indices.end(); ++i) {
      auto ii = i; ++ii;
      stride.push_back(std::accumulate(ii, indices.end(), 1ul, prod));
    }

    std::vector<size_t> extent(rank);
    auto e = extent.rbegin();
    for (int i = 0; i != rank; ++i)
      *e++ = input->extent(i);

    std::vector<size_t> stride_target;
    for (auto i = extent.begin(); i != extent.end(); ++i) {
      auto ii = i; ++ii;
      stride_target.push_back(std::accumulate(ii, extent.end(), 1ul, std::multiplies<size_t>()));
    }

    typename TiledArray::Array<DataType,N>::value_type tile(range);
    const size_t buffersize = std::accumulate(indices.begin(), indices.end(), 1ul, prod);
    const size_t backsize = indices.back().size();
    for (size_t n = 0; n != buffersize; n += backsize) {
      size_t offset = 0lu;
      size_t tmp = n;
      for (int i = 0; i != rank; ++i) {
        offset += (tmp / stride[i] + indices[i].offset() - inpoffsets[i]) * stride_target[i];
        tmp = n % stride[i];
      }
      std::copy_n(input->data()+offset, backsize, &(tile[n]));
    }
    *it = tile;
  }
}


}
}

#endif

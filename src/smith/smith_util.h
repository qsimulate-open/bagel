//
// BAGEL - Parallel electron correlation program.
// Filename: smith_util.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_SMITH_SMITH_UTIL_H
#define __SRC_SMITH_SMITH_UTIL_H

#include <src/smith/tensor.h>
#include <src/util/kramers.h>

namespace bagel {
namespace SMITH {

template<int N, typename DataType>
static void fill_block(std::shared_ptr<Tensor_<DataType>> target, std::shared_ptr<const btas::TensorN<DataType,N>> input,
                       const std::vector<int>& inpoffsets, const std::vector<IndexRange>& ranges_rev) {
  assert(input->range().ordinal().contiguous());
  assert(target->rank() == input->range().rank() && target->rank() > 0);
  const int rank = target->rank();
  const std::vector<IndexRange> ranges(ranges_rev.rbegin(), ranges_rev.rend());

  auto prod = [](const size_t n, const Index& i) { return n*i.size(); };

  LoopGenerator gen(ranges);
  std::vector<std::vector<Index>> loop = gen.block_loop();
  for (auto& indices : loop) {
    assert(indices.size() == rank);

    const size_t buffersize = std::accumulate(indices.begin(), indices.end(), 1ul, prod);
    std::unique_ptr<DataType[]> buffer(new DataType[buffersize]);
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

    const size_t backsize = indices.back().size();
    for (size_t n = 0; n != buffersize; n += backsize) {
      size_t offset = 0lu;
      size_t tmp = n;
      for (int i = 0; i != rank; ++i) {
        offset += (tmp / stride[i] + indices[i].offset() - inpoffsets[i]) * stride_target[i];
        tmp = n % stride[i];
      }
      std::copy_n(input->data()+offset, backsize, buffer.get()+n);
    }

    target->put_block(buffer, std::vector<Index>(indices.rbegin(), indices.rend()));
  }
}


template<int N, typename DataType>
static void fill_block(std::shared_ptr<TATensor<DataType,N>> target, std::shared_ptr<const btas::TensorN<DataType,N>> input,
                       const std::vector<int>& inpoffsets) {
  assert(input->range().ordinal().contiguous());
  assert(N == input->range().rank());
  const int rank = N;

  auto prod = [](const size_t n, const Index& i) { return n*i.size(); };

  // first index corresponds to the fastest in the data
  const std::vector<IndexRange> ranges_rev = target->indexrange();
  // last index corresponds to the fastest
  const std::vector<IndexRange> ranges(ranges_rev.rbegin(), ranges_rev.rend());

  // The last element corresponds to the fastest index
  std::vector<std::map<size_t,Index>> keymap;
  for (auto it = ranges.begin(); it != ranges.end(); ++it) {
    std::map<size_t,Index> key;
    for (auto& j : *it)
      key.emplace(j.offset()-it->front().offset(), j);
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


template<int N, typename DataType, class T> // T is supposed to be derived from btas::Tensor
static void fill_block(std::shared_ptr<Tensor_<DataType>> target, std::shared_ptr<const Kramers<N,T>> input,
                       const std::vector<int>& inpoffsets, const std::vector<IndexRange>& ranges_rev) {
  const int rank = target->rank();
  const std::vector<IndexRange> ranges(ranges_rev.rbegin(), ranges_rev.rend());

  auto prod = [](const size_t n, const Index& i) { return n*i.size(); };

  LoopGenerator gen(ranges);
  std::vector<std::vector<Index>> loop = gen.block_loop();
  for (auto& indices : loop) {
    assert(indices.size() == rank);

    std::bitset<N> bit;
    for (int i = 0; i != N; ++i)
      bit[i] = indices[i].kramers() ? 1 : 0; // indices is reversed, so this is correct

    if (input->exist(bit)) {
      const size_t buffersize = std::accumulate(indices.begin(), indices.end(), 1ul, prod);
      std::vector<size_t> stride;
      for (auto i = indices.begin(); i != indices.end(); ++i) {
        auto ii = i; ++ii;
        stride.push_back(std::accumulate(ii, indices.end(), 1ul, prod));
      }

      std::unique_ptr<DataType[]> buffer(new DataType[buffersize]);

      // in principle there is repetition (especially when active orbitals are separated into small blocks)
      auto block = input->at(bit);
      assert(block->range().ordinal().contiguous());
      assert(target->rank() == block->range().rank() && target->rank() > 0);

      std::vector<size_t> extent(rank);
      auto e = extent.rbegin();
      for (int i = 0; i != rank; ++i)
        *e++ = block->extent(i);

      std::vector<size_t> stride_target;
      for (auto i = extent.begin(); i != extent.end(); ++i) {
        auto ii = i; ++ii;
        stride_target.push_back(std::accumulate(ii, extent.end(), 1ul, std::multiplies<size_t>()));
      }

      const size_t backsize = indices.back().size();
      for (size_t n = 0; n != buffersize; n += backsize) {
        size_t offset = 0lu;
        size_t tmp = n;
        for (int i = 0; i != rank; ++i) {
          offset += (tmp / stride[i] + indices[i].kramers_offset() - inpoffsets[i]) * stride_target[i];
          tmp = n % stride[i];
        }
        std::copy_n(block->data()+offset, backsize, buffer.get()+n);
      }
      target->put_block(buffer, std::vector<Index>(indices.rbegin(), indices.rend()));
    }
  }
  target->set_perm(input->perm());
}

}
}

#endif

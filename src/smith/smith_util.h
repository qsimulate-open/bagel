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
static std::shared_ptr<Tensor_<DataType>>
  fill_block(std::shared_ptr<const btas::TensorN<DataType,N>> input, const std::vector<int>& inpoffsets_rev, const std::vector<IndexRange>& ranges_rev) {

  auto target = std::make_shared<Tensor_<DataType>>(ranges_rev);
  target->allocate();

  assert(input->range().ordinal().contiguous());
  assert(target->rank() == input->range().rank() && target->rank() > 0);
  const int rank = target->rank();
  const std::vector<IndexRange> ranges(ranges_rev.rbegin(), ranges_rev.rend());
  const std::vector<int> inpoffsets(inpoffsets_rev.rbegin(), inpoffsets_rev.rend());

  auto prod = [](const size_t n, const Index& i) { return n*i.size(); };

  std::vector<std::vector<Index>> loop = LoopGenerator::gen(ranges);
  for (auto& indices : loop) {
    assert(indices.size() == rank);
    if (!target->is_local(std::vector<Index>(indices.rbegin(), indices.rend()))) continue;

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
  mpi__->barrier();
  return target;
}


template<int N, typename DataType, class T> // T is supposed to be derived from btas::Tensor
static std::shared_ptr<Tensor_<DataType>>
  fill_block(std::shared_ptr<const Kramers<N,T>> input, const std::vector<int>& inpoffsets_rev, const std::vector<IndexRange>& ranges_rev) {

  const int rank = N;
  const std::vector<IndexRange> ranges(ranges_rev.rbegin(), ranges_rev.rend());
  const std::vector<int> inpoffsets(inpoffsets_rev.rbegin(), inpoffsets_rev.rend());
  const std::list<std::vector<bool>> kramers = input->listkramers();

  std::vector<std::vector<Index>> loop = LoopGenerator::gen(ranges);

  std::unordered_set<size_t> sparse;
  for (auto& indices : loop) {
    std::vector<bool> tmp;
    for (auto i = indices.rbegin(); i != indices.rend(); ++i)
      tmp.push_back(i->kramers());
    if (std::find(kramers.begin(), kramers.end(), tmp) != kramers.end())
    sparse.insert(generate_hash_key(std::vector<Index>(indices.rbegin(), indices.rend())));
  }

  auto target = std::make_shared<Tensor_<DataType>>(ranges_rev, true, sparse, true);
  target->set_stored_sectors(kramers);

  auto prod = [](const size_t n, const Index& i) { return n*i.size(); };

  for (auto& indices : loop) {
    assert(indices.size() == rank);

    std::bitset<N> bit;
    for (int i = 0; i != N; ++i)
      bit[i] = indices[i].kramers() ? 1 : 0; // indices is reversed, so this is correct

    if (input->exist(bit)) {
      if (!target->is_local(std::vector<Index>(indices.rbegin(), indices.rend()))) continue;

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
  mpi__->barrier();
  return target;
}

}
}

#endif

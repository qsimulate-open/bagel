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

}
}

#endif

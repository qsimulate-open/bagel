//
// BAGEL - Parallel electron correlation program.
// Filename: tatensor.h
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

#ifndef __SRC_SMITH_TATENSOR_H
#define __SRC_SMITH_TATENSOR_H

#include <tiledarray.h>
#include <src/smith/indexrange.h>

namespace bagel {
namespace SMITH {

// wrapper class for TiledArray Tensor
template<typename DataType, int N>
class TATensor {
  protected:
    std::shared_ptr<TiledArray::Array<DataType,N>> data_;
    std::vector<IndexRange> range_;

  public:
    TATensor(const std::vector<IndexRange>& r) : range_(r) {
      assert(r.size() == N);
      std::vector<TiledArray::TiledRange1> ranges;
      for (auto it = range_.rbegin(); it != range_.rend(); ++it) {
        std::vector<size_t> tile_boundaries;
        for (auto& j : *it)
          tile_boundaries.push_back(j.offset());
        tile_boundaries.push_back(it->back().offset()+it->back().size());
        ranges.emplace_back(tile_boundaries.begin(), tile_boundaries.end());
      }
      TiledArray::TiledRange trange(ranges.begin(), ranges.end());
      data_ = std::make_shared<TiledArray::Array<DataType, N>>(madness::World::get_default(), trange);
    }

    TATensor(const TATensor<DataType,N>& o) : data_(std::make_shared<TATensor<DataType,N>>(o.data_)), range_(o.range_) {
    }

    TATensor(TATensor<DataType,N>&& o) : data_(std::make_shared<TATensor<DataType,N>>(std::move(o.data_))), range_(o.range_) {
    }

    std::shared_ptr<TiledArray::Array<DataType,N>> data() { return data_; }
    std::shared_ptr<const TiledArray::Array<DataType,N>> data() const { return data_; }

    std::vector<IndexRange> indexrange() const { return range_; }
};

}}

#endif

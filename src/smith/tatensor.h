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

//TODO do something for symmetry
//  std::map<std::vector<int>, std::pair<double,bool>> perm_;

  public:
    TATensor(const std::vector<IndexRange>& r, /*toggle for symmetry*/const bool dummy = false) : range_(r) {
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

    auto local(const std::vector<Index>& index) -> decltype(std::make_pair(true,data_->begin())) {
      assert(index.size() == N);
      bool out = false;
      auto it = data_->begin();
      for ( ; it != data_->end(); ++it) {
        const TiledArray::Range range = data_->trange().make_tile_range(it.ordinal());
        auto lo = range.lobound();
        assert(lo.size() == N);
        bool found = true;
        auto j = index.rbegin();
        for (auto& i : lo) {
          found &= i == j->offset();
          ++j;
        }
        out = found;
        if (found) break;
      }
      return std::make_pair(out, it);
    }

    std::shared_ptr<TiledArray::Array<DataType,N>> data() { return data_; }
    std::shared_ptr<const TiledArray::Array<DataType,N>> data() const { return data_; }

    std::vector<IndexRange> indexrange() const { return range_; }
};

}}

#endif

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: overlap_point.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_PROP_OVERLAP_POINT_H
#define __SRC_PROP_OVERLAP_POINT_H

#include <src/wfn/geometry.h>
#include <src/integral/compos/point_complexoverlapbatch.h>
#include <src/integral/os/point_overlapbatch.h>

namespace bagel {

template <typename MatType, typename IntType,
          class Enable = typename std::enable_if<((std::is_same<MatType, Matrix>::value && std::is_same<IntType, Point_OverlapBatch>::value) ||
                                                  (std::is_same<MatType, ZMatrix>::value && std::is_same<IntType, Point_ComplexOverlapBatch>::value))>::type>
class Overlap_Point_ {
  private:
    using DataType = typename MatType::value_type;

  protected:
    std::shared_ptr<const Geometry> geom_;
    std::array<double,3> location_;

  public:
    Overlap_Point_(std::shared_ptr<const Geometry> g, std::array<double,3> loc);

    std::shared_ptr<MatType> compute() const;
};

using Overlap_Point = Overlap_Point_<Matrix, Point_OverlapBatch>;
using Overlap_Point_London = Overlap_Point_<ZMatrix, Point_ComplexOverlapBatch>;

}

extern template class bagel::Overlap_Point_<bagel::Matrix, bagel::Point_OverlapBatch>;
extern template class bagel::Overlap_Point_<bagel::ZMatrix, bagel::Point_ComplexOverlapBatch>;

#endif

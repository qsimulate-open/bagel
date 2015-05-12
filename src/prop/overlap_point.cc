//
// BAGEL - Parallel electron correlation program.
// Filename: overlap_point.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/prop/overlap_point.h>

using namespace std;
using namespace bagel;

template <typename MatType, typename IntType, class Enable>
Overlap_Point_<MatType, IntType, Enable>::Overlap_Point_(shared_ptr<const Geometry> g, const array<double,3> loc) : geom_(g), location_(loc) { }


template <typename MatType, typename IntType, class Enable>
shared_ptr<MatType> Overlap_Point_<MatType, IntType, Enable>::compute() const {

  const int nbasis = geom_->nbasis();
  auto out = make_shared<MatType>(nbasis, nbasis);

  // TODO perhaps we could reduce operation by a factor of 2
  auto o0 = geom_->offsets().begin();
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++o0) {
    auto o1 = geom_->offsets().begin();
    for (auto a1 = geom_->atoms().begin(); a1 != geom_->atoms().end(); ++a1, ++o1) {

      auto offset0 = o0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++offset0) {
        auto offset1 = o1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++offset1) {

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          IntType batch(input, geom_->magnetic_field(), location_);
          batch.compute();

          const DataType* dat0 = batch.data();
          for (int i = *offset0; i != *offset0 + (*b0)->nbasis(); ++i) {
            for (int j = *offset1; j != *offset1 + (*b1)->nbasis(); ++j, ++dat0) {
              out->element(j,i) = *dat0;
            }
          }

        }
      }
    }
  }

  return out;
}


template class Overlap_Point_<bagel::Matrix, Point_OverlapBatch>;
template class Overlap_Point_<bagel::ZMatrix, Point_ComplexOverlapBatch>;


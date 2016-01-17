//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: point_overlapbatch.h
// Copyright (C) 2015 Toru Shiozaki
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

// Analogous to OverlapBatch, but instead of integrating over all space,
// we just get the value at location_.  This is used to visualize molecular orbitals.

#ifndef __SRC_INTEGRAL_COMPOS_POINT_OVERLAPBATCH_H
#define __SRC_INTEGRAL_COMPOS_POINT_OVERLAPBATCH_H

#include <src/integral/os/osintegral.h>

namespace bagel {

class Point_OverlapBatch : public OSIntegral<double,Int_t::Standard> {
  protected:
    std::array<double,3> location_;

    void perform_VRR(double*) override;

    int nblocks() const override { return 1; }
    int nrank() const override { return 0; }

    bool diagonal_;

  public:
    // TODO Can dummy argument be avoided?
    Point_OverlapBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, const std::array<double,3> dummy, const std::array<double,3> _loc,
                         std::shared_ptr<StackMem> stack = nullptr)
     : OSIntegral<double,Int_t::Standard>(basis, stack), location_(_loc) { common_init();
       if (*basisinfo_[0] == *basisinfo_[1])
         diagonal_ = true;
       else
         diagonal_ = false;
      }

    void compute() override;
};

}

#endif

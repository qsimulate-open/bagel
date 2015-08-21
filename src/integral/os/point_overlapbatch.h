//
// BAGEL - Parallel electron correlation program.
// Filename: point_overlapbatch.h
// Copyright (C) 2015 Toru Shiozaki
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

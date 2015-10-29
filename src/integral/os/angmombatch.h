//
// BAGEL - Parallel electron correlation program.
// Filename: angmombatch.h
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

// Orbital angular momentum integrals, for real Gaussian basis functions, around some center mcoord_
// The result is pure imaginary, so we actually return <A| -iL |B>; multiply by i for the correct matrix elements

#ifndef __SRC_INTEGRAL_OS_ANGMOMBATCH_H
#define __SRC_INTEGRAL_OS_ANGMOMBATCH_H

#include <src/integral/os/osintegral.h>
#include <src/molecule/molecule.h>

namespace bagel {

class AngMomBatch : public OSIntegral<double,Int_t::Standard> {
  protected:
    std::array<double,3> mcoord_;

    void perform_VRR(double*) override;

    int nblocks() const override { return 3; }
    int nrank() const override { return 0; }

  public:
    AngMomBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, const std::array<double,3> _mcoord)
      : OSIntegral<double,Int_t::Standard>(basis), mcoord_(_mcoord) { common_init(); }

    constexpr static int Nblocks() { return 3; }

    void compute() override;
};

}

#endif

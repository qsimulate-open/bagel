//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: angmombatch.h
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

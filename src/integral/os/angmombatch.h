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

// Orbital angular momentum integrals, for GIAO basis functions, around some center mcoord_

#ifndef __SRC_INTEGRAL_OS_ANGMOMBATCH_H
#define __SRC_INTEGRAL_OS_ANGMOMBATCH_H

#include <src/integral/os/osintegral.h>

namespace bagel {

class AngMomBatch : public OSIntegral<std::complex<double>,Int_t::London> {
  protected:
    std::array<double,3> mcoord_;

    void perform_VRR(std::complex<double>*) override;

    int nblocks() const override { return 3; }
    int nrank() const override { return 0; }

  public:
    AngMomBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, const std::array<double,3> _mcoord)
      : OSIntegral<std::complex<double>,Int_t::London>(basis), mcoord_(_mcoord) { common_init(); }

    void compute() override;
};

}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: jexpansion.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_PERIODIC_JEXPANSION_H
#define __SRC_PERIODIC_JEXPANSION_H

#include <src/periodic/multipolebatch.h>
#include <src/periodic/localexpansion.h>

namespace bagel {

// compute (rs|tu) using multipole expansion if rs and tu are well-separated
class JExpansion {
  protected:
    std::array<std::shared_ptr<const Shell>, 4> basisinfo_;
    int lmax_, ws_;

    std::array<double, 3> centre0_, centre1_, r12_;
    double extent0_, extent1_;
    int num_multipoles_;

    std::array<double, 3> distribution_centre(std::array<std::shared_ptr<const Shell>, 2> shells);
    double distribution_extent(std::array<std::shared_ptr<const Shell>, 2> shells, const double thresh = PRIM_SCREEN_THRESH);

  public:
    JExpansion(const std::array<std::shared_ptr<const Shell>,4>& shells, const int lmax = ANG_HRR_END, const int ws = 2);
    ~JExpansion() { }

    void init();
    bool is_well_separated();
    std::vector<std::shared_ptr<const ZMatrix>> compute(std::shared_ptr<const Matrix> density);
};

}

#endif

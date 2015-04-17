//
// BAGEL - Parallel electron correlation program.
// Filename: localexpansion.h
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


#ifndef __SRC_PERIODIC_LOCALEXPANSION_H
#define __SRC_PERIODIC_LOCALEXPANSION_H

#include <vector>
#include <src/util/constants.h>
#include <src/util/math/legendre.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

class LocalExpansion {
  protected:
    std::array<double, 3> centre_;
    std::vector<std::shared_ptr<const ZMatrix>> moments_;
    int lmax_;

    int nbasis0_, nbasis1_;
    int num_multipoles_;

  public:
    LocalExpansion(const std::array<double, 3> centre, std::vector<std::shared_ptr<const ZMatrix>> moments,
                   const int lmax = ANG_HRR_END);
    ~LocalExpansion() { }

    int lmax() const { return lmax_; }
    std::array<double, 3> centre() const { return centre_; }
    double centre(const int i) const { return centre_[i]; }

    std::vector<std::shared_ptr<const ZMatrix>> moments() const { return moments_; }
    std::shared_ptr<const ZMatrix> moment(const int i) const { return moments_[i]; }
    std::shared_ptr<const ZMatrix> moment(const int l, const int m) const { return moments_[l * l + l + m]; }

    std::vector<std::shared_ptr<const ZMatrix>> compute_local_moments();
    std::vector<std::shared_ptr<const ZMatrix>> compute_shifted_moments();
    std::vector<std::shared_ptr<const ZMatrix>> compute_shifted_local_expansions();
};

}

#endif

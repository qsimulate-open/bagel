//
// BAGEL - Parallel electron correlation program.
// Filename: angularbatch.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_INTEGRAL_ANGULARBATCH_H
#define __SRC_INTEGRAL_ANGULARBATCH_H

#include <src/math/bessel.h>
#include <src/util/constants.h>
#include <src/integral/ecp/radial.h>
#include <src/integral/ecp/sphharmonics.h>
#include <src/integral/ecp/wigner3j.h>
#include <src/integral/ecp/ecpbatch.h>

namespace bagel {

class AngularBatch: public ECPBatch {
  protected:

    double integrate3SHs(std::array<std::pair<int, int>, 3> lm) const;
    double integrate3USP(std::array<int, 3> xyz_exponents) const;
    double integrate2SH1USP(const std::pair<int, int> lm1, const std::pair<int, int> lm2, const std::array<int, 3> ijk) const;
    double project_one_gaussian(std::array<double, 3> posA, std::array<int, 3> lxyz, const double expA,
                                std::array<double, 3> posB, std::array<int, 2> lm, const double r) const;

    double* projdata_;

  public:
    AngularBatch(const std::shared_ptr<const Shell_ECP>& _ecp_info, const std::array<std::shared_ptr<const Shell>,2>& _info,
            const std::shared_ptr<const Molecule> mol)
     : ECPBatch(_ecp_info, _info, mol) {
//   const double integral_thresh = PRIM_SCREEN_THRESH;
      /// call some functions to compute integrals
    }

    ~AngularBatch() {}

};

}

#endif

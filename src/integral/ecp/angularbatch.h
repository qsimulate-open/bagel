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


#ifndef __SRC_INTEGRAL_ECP_ANGULARBATCH_H
#define __SRC_INTEGRAL_ECP_ANGULARBATCH_H

#include <src/util/constants.h>
#include <src/molecule/atom.h>
#include <src/integral/ecp/sphharmonics.h>

namespace bagel {

class AngularBatch {
  protected:

    std::array<std::shared_ptr<const Shell>,2> basisinfo_;
    std::shared_ptr<const ECP> ecp_;
    double exp0_, exp1_;
    std::array<int, 3> ang0_, ang1_;

    double integrate3SHs(std::array<std::pair<int, int>, 3> lm) const;
    double integrate3USP(std::array<int, 3> xyz_exponents) const;
    double integrate2SH1USP(const std::pair<int, int> lm1, const std::pair<int, int> lm2, const std::array<int, 3> ijk) const;
    double project_one_centre(std::array<double, 3> posA, const std::array<int, 3> lxyz, const double expA,
                              std::array<double, 3> posB, const std::array<int, 2> lm, const double r);
    double project_many_centres(const double expA, const double expC, const double r);

  public:
    AngularBatch(const std::shared_ptr<const ECP> _ecp, const std::array<std::shared_ptr<const Shell>,2>& _info,
                 const double expA, const double expC,
                 const std::array<int, 3> angA, const std::array<int, 3> angC)
     : basisinfo_(_info), ecp_(_ecp), exp0_(expA), exp1_(expC), ang0_(angA), ang1_(angC) {}

    ~AngularBatch() {}
    double compute(const double r);

    void print() const;
    void print_one_centre(std::array<double, 3> posA, const std::array<int, 3> lxyz, const double expA,
                          std::array<double, 3> posB, const std::array<int, 2> lm, const double r) const;

};

}

#endif

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gsmallnaibatch.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_INTEGRAL_RYS_GSMALLNAIBATCH_H
#define __SRC_INTEGRAL_RYS_GSMALLNAIBATCH_H

#include <src/util/math/xyzfile.h>
#include <src/integral/smallints1e.h>

namespace bagel {

// TODO make a base class SmallNAIBatch_base and remove redundancy
class GSmallNAIBatch {
  protected:
    std::vector<std::shared_ptr<Matrix>> data_;

    const std::shared_ptr<const Molecule> mol_;
    const std::array<std::shared_ptr<const Shell>,2> shells_;

    const size_t size_block_;

    // specific to gradients
    std::tuple<int,int> iatom_;

  public:
    GSmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Molecule> mol, const std::tuple<int,int> i);

    // computes derivative NAI over Cartesian (i.e., RI) basis functions.
    void compute();

    // multiplies S^-1P to density and contracts with cartesian integrals.
    std::shared_ptr<GradFile> compute_gradient(std::array<std::shared_ptr<const Matrix>,6> d) const;

    size_t size_block() const { return size_block_; }

};

}

#endif

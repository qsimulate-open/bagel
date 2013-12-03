//
// BAGEL - Parallel electron correlation program.
// Filename: gsmallnaibatch.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_INTEGRAL_RYS_GSMALLNAIBATCH_H
#define __SRC_INTEGRAL_RYS_GSMALLNAIBATCH_H

#include <src/math/xyzfile.h>
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

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
// the Free Software Foundation; either version 2, or (at your option)
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

#include <src/grad/gradfile.h>
#include <src/integral/rys/smallnaibatch.h>

namespace bagel {

class GSmallNAIBatch : public SmallNAIBatch {
  protected:
    std::tuple<int,int> iatom_;

    std::shared_ptr<StackMem> stack_;
    bool allocated_here_;

  public:
    GSmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Molecule> mol, const std::tuple<int,int> i,
                   std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>())
    : SmallNAIBatch(info, mol), iatom_(i) {
      stack_ = stack ? stack : resources__->get();
      allocated_here_ = stack != stack_;
    }

    ~GSmallNAIBatch() {
      if (allocated_here_)
        resources__->release(stack_);
    }

    // computes derivative NAI over Cartesian (i.e., RI) basis functions.
    void compute() override {}

    // multiplies S^-1P to density and contracts with cartesian integrals.
    std::shared_ptr<GradFile> compute_gradient(std::array<std::shared_ptr<const Matrix>,6> d) const { return std::make_shared<GradFile>(mol_->natom()); }

};

}

#endif

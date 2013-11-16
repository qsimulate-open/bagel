//
// BAGEL - Parallel electron correlation program.
// Filename: sonaibatch.h
// Copyright (C) 2009 Toru Shiozaki
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

#ifndef __SRC_INTEGRAL_RYS_SONAIBATCH_H
#define __SRC_INTEGRAL_RYS_SONAIBATCH_H

#include <src/integral/rys/sonaibatch_base.h>

namespace bagel {

class SONAIBatch : public SONAIBatch_base {

  protected:

  public:

    SONAIBatch(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>())
      :  SONAIBatch_base(_info, mol, 0, stack, 0, 0.0) {};

    SONAIBatch(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, const int L, const double A = 0.0)
      :  SONAIBatch_base(_info, mol, 0, std::shared_ptr<StackMem>(), L, A) {};
     ~SONAIBatch() {};

    /// compute a batch of integrals
    void compute() override;

    constexpr static int Nblocks() { return 1; }
};

}

#endif


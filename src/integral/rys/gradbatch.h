//
// BAGEL - Parallel electron correlation program.
// Filename: gradbatch.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_RYSINT_GRADBATCH_H
#define __SRC_RYSINT_GRADBATCH_H

// compute analytic nuclear gradients
#include <stddef.h>
#include <memory>
#include <src/integral/rys/eribatch_base.h>

namespace bagel {

class GradBatch : public ERIBatch_base {
  protected:
    // if we only compute three-center integrals, we want to use this info
    // to reduce the number of differentiation
    int centers_;

    void perform_VRR();

    void set_exponents();
    std::unique_ptr<double[]> exponents_;

  public:
    GradBatch(const std::array<std::shared_ptr<const Shell>,4>& shells, const double max_density, const double dummy = 0.0, const bool dum = true,
              std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>());

    void compute() override;

    double* data(const int i) override { assert(i < 12); return data_ + i*size_block_; }

};

}

#endif


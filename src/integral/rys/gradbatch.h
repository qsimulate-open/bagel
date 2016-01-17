//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gradbatch.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_RYSINT_GRADBATCH_H
#define __SRC_RYSINT_GRADBATCH_H

// compute analytic nuclear gradients
#include <src/integral/rys/eribatch_base.h>

namespace bagel {

class GradBatch : public ERIBatch_base {
  protected:
    // if we only compute three-center integrals, we want to use this info
    // to reduce the number of differentiation
    int centers_;

    void perform_VRR();
    void root_weight(const int ps) override;

    void set_exponents();
    std::unique_ptr<double[]> exponents_;

  public:
    GradBatch(const std::array<std::shared_ptr<const Shell>,4>& shells, const double max_density, const double dummy = 0.0, const bool dum = true,
              std::shared_ptr<StackMem> stack = nullptr);

    void compute() override;

    double* data(const int i) override { assert(i < 12); return data_ + i*size_block_; }

};

}

#endif


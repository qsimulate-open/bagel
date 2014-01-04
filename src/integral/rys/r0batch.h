//
// BAGEL - Parallel electron correlation program.
// Filename: r0batch.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_RYSINT_R0BATCH_H
#define __SRC_RYSINT_R0BATCH_H

#include <src/integral/rys/rnbatch.h>

namespace bagel {

class R0Batch: public RnBatch {
  protected:
    void compute_ssss(const double) override;
    double scale_root(const double root, const double p, const double zeta) override { return -zeta / (p - zeta); } //zeta = -zeta for testing
    double scale_weight(const double weight, const double coef) override { return coef; }

  public:
    R0Batch(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol,
            std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>())
      : RnBatch (_info, mol, stack) {
      rank_ = 1;
      const double integral_thresh = PRIM_SCREEN_THRESH;
      compute_ssss(integral_thresh);
    }


    R0Batch(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol,
            const int L, const double A = 0.0)
      : RnBatch (_info, mol, L, A) {
      rank_ = 1;
      const double integral_thresh = PRIM_SCREEN_THRESH;
      compute_ssss(integral_thresh);
    }


    ~R0Batch() {}

};

}

#endif

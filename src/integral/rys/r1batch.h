//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: r1batch.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_RYSINT_R1BATCH_H
#define __SRC_RYSINT_R1BATCH_H

#include <src/integral/rys/rnbatch.h>

namespace bagel {

class R1Batch: public RnBatch {
  protected:

    void root_weight(const int ps) override;
    void compute_ssss(const double) override;
    double scale_root(const double root, const double p, const double zeta) override { return (p * root + zeta)/(p + zeta); }

  public:
    R1Batch(const std::array<std::shared_ptr<const Shell>,2>& _info,
            const std::shared_ptr<const Molecule> mol, std::shared_ptr<StackMem> stack = nullptr)
      : RnBatch (_info, mol, stack) {
      const double integral_thresh = PRIM_SCREEN_THRESH;
      compute_ssss(integral_thresh);
      root_weight(primsize_*natom_*max_rterms_);
    }

    ~R1Batch() {}

};

}

#endif

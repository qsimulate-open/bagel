//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: breitbatch.h
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

#ifndef __SRC_INTEGRAL_RYS_BREITBATCH_H
#define __SRC_INTEGRAL_RYS_BREITBATCH_H

#include <src/integral/rys/eribatch_base.h>

namespace bagel {

class BreitBatch_base : public ERIBatch_base {
  protected:
    virtual void perform_VRR() = 0;
    virtual void perform_VRR1() = 0;
  public:
    BreitBatch_base(const std::array<std::shared_ptr<const Shell>,4>& a, const double max_density, const double dummy, const bool dum, const int N)
      :  ERIBatch_base(a, 0, N) {

      const double integral_thresh = (max_density != 0.0) ? (PRIM_SCREEN_THRESH / max_density) : 0.0;
      compute_ssss(integral_thresh);

      root_weight(this->primsize_);
    }

    /// compute a batch of integrals
    void compute() override;

    void root_weight(const int ps) override;
    double* data(const int i) override { return data_ + i*size_block_; }
    constexpr static int nblocks() { return 6; }
};


class BreitBatch : public BreitBatch_base {
  protected:
    void perform_VRR() override;
    void perform_VRR1() override;

  public:
    // dummy will never be used.
    BreitBatch(const std::array<std::shared_ptr<const Shell>,4>& a, const double max_density, const double dummy = 0.0, const bool dum = true)
     :  BreitBatch_base(a, max_density, dummy, dum, 1) { }

};


class Spin2Batch : public BreitBatch_base {
  protected:
    void perform_VRR() override;
    void perform_VRR1() override;

  public:
    // dummy will never be used.
    Spin2Batch(const std::array<std::shared_ptr<const Shell>,4>& a, const double max_density, const double dummy = 0.0, const bool dum = true)
     :  BreitBatch_base(a, max_density, dummy, dum, 2) { }

};

}

#endif


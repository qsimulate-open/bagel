//
// BAGEL - Parallel electron correlation program.
// Filename: spin2batch.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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

#ifndef __SRC_INTEGRAL_RYS_SPIN2BATCH_H
#define __SRC_INTEGRAL_RYS_SPIN2BATCH_H

#include <src/integral/rys/breitbatch.h>

namespace bagel {

class Spin2Batch : public BreitBatch {

  protected:
    void perform_VRR();
    void perform_VRR1();

  public:

    // dummy will never used.
    Spin2Batch(const std::array<std::shared_ptr<const Shell>,4>&, const double max_density, const double dummy = 0.0, const bool dum = true);

    /// compute a batch of integrals
    void compute() override;

    double* data(const int i) override { return data_ + i*size_block_; }
    constexpr static int nblocks() { return 6; }

};

}

#endif


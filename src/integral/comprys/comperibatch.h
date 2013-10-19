#if 0
//
// BAGEL - Parallel electron correlation program.
// Filename: comperibatch.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Ryan Reynolds <rreynoldschem@u.northwestern.edu>
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

#ifndef __SRC_INTEGRAL_COMPRYS_COMPERIBATCH_H
#define __SRC_INTEGRAL_COMPRYS_COMPERIBATCH_H

#include <src/integral/rys/eribatch_base.h>
#include <complex>

namespace bagel {

class CompERIBatch : public ERIBatch_Base<complex<double>> {

  protected:
    void perform_VRR();
    void perform_VRR1();
    void perform_VRR2();
    void perform_VRR3();

  public:

    // dummy will never be used.
    CompERIBatch(const std::array<std::shared_ptr<const Shell>,4>&, const double max_density, const double dummy = 0.0, const bool dum = true,
             std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>());

    /// compute a batch of integrals
    void compute() override;

    constexpr static int Nblocks() { return 1; }
};

}

#endif
#endif

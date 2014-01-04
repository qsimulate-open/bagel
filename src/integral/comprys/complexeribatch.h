//
// BAGEL - Parallel electron correlation program.
// Filename: complexeribatch.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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

#ifndef __SRC_INTEGRAL_COMPRYS_COMPLEXERIBATCH_H
#define __SRC_INTEGRAL_COMPRYS_COMPLEXERIBATCH_H

#include <src/integral/rys/eribatch_base.h>
#include <complex>

namespace bagel {

class ComplexERIBatch : public ERIBatch_Base<std::complex<double>,Int_t::London> {

  protected:
    void perform_VRR();
    void root_weight(const int ps) override;
    std::complex<double> get_PQ(const double coord1, const double coord2, const double exp1, const double exp2, const double one12, const int center1, const int dim, const bool swap) override;

  public:

    // dummy will never be used.
    ComplexERIBatch(const std::array<std::shared_ptr<const Shell>,4>&, const double max_density, const std::complex<double> dummy = std::complex<double>(0.0,0.0),
             const bool dum = true, std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>());

    /// compute a batch of integrals
    void compute() override;

    constexpr static int Nblocks() { return 1; }

};

}

#endif

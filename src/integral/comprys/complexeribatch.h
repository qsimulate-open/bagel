//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexeribatch.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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
             const bool dum = true, std::shared_ptr<StackMem> stack = nullptr);

    /// compute a batch of integrals
    void compute() override;

    constexpr static int Nblocks() { return 1; }

};

}

#endif

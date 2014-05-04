//
// BAGEL - Parallel electron correlation program.
// Filename: complexnaibatch.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#ifndef __SRC_INTEGRAL_COMPRYS_COMPLEXNAIBATCH_H
#define __SRC_INTEGRAL_COMPRYS_COMPLEXNAIBATCH_H

#include <src/integral/sortlist.h>
#include <src/integral/carsphlist.h>
#include <src/integral/hrrlist.h>
#include <src/integral/rys/coulombbatch_base.h>

namespace bagel {

class ComplexNAIBatch : public CoulombBatch_Base<std::complex<double>,Int_t::London> {

  protected:

    void root_weight(const int ps) override;
    std::complex<double> get_PQ(const double coord1, const double coord2, const double exp1, const double exp2, const double one12, const int center1, const int dim, const bool swap) override;

  public:

    void compute() override;

    ComplexNAIBatch(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, std::shared_ptr<StackMem> stack = nullptr)
      :  CoulombBatch_Base<std::complex<double>,Int_t::London>(_info, mol, 0, stack, 0, 0.0) {
      const double integral_thresh = PRIM_SCREEN_THRESH;
      compute_ssss(integral_thresh);
      root_weight(primsize_*natom_);
    }

    ComplexNAIBatch(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, const int L, const double A = 0.0)
      :  CoulombBatch_Base<std::complex<double>,Int_t::London>(_info, mol, 0, nullptr, L, A) {
      const double integral_thresh = PRIM_SCREEN_THRESH;
      compute_ssss(integral_thresh);
      root_weight(primsize_*natom_);
    }

    ~ComplexNAIBatch() {};

    constexpr static int Nblocks() { return 1; }

};

}

#endif


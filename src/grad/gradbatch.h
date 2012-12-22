//
// BAGEL - Parallel electron correlation program.
// Filename: gradbatch.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#ifndef __SRC_RYSINT_GRAD_H
#define __SRC_RYSINT_GRAD_H

// compute analytic nuclear gradients
#include <stddef.h>
#include <memory>
#include <src/rysint/eribatch_base.h>

namespace bagel {

class GradBatch : public ERIBatch_base {
  protected:
    // if we only compute three-center integrals, we want to use this info
    // to reduce the number of differentiation
    int centers_;

    // in the gradient computation, we don't assemble until we operate differential operators.
    // we don't use translational invariance so far.
    void perform_VRR1();
    void perform_VRR2();
    void perform_VRR3();
    void perform_VRR4();
    void perform_VRR5();
    void perform_VRR6();
    void perform_VRR7();
    void perform_VRR8();
    void perform_VRR9();
    void perform_VRR10();
    void perform_VRR11();
    void perform_VRR12();
    void perform_VRR13();

    void set_exponents();
    std::unique_ptr<double[]> exponents_;

    template<int rank>
    size_t m(const int& i, const int& a, const int& b, const int& c, const int& d) const {
      const int la = basisinfo_[0]->angular_number()+2;
      const int lb = basisinfo_[1]->angular_number()+2;
      const int lc = basisinfo_[2]->angular_number()+2;
      const int ld = basisinfo_[3]->angular_number()+2;
      return i+rank*(a+la*(b+lb*(c+lc*d)));
    }

  public:
    GradBatch(const std::array<std::shared_ptr<const Shell>,4>& shells, const double max_density, const double dummy = 0.0, const bool dum = true);
    ~GradBatch();

    void compute() override;

};

}

#endif


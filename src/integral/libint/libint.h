//
// BAGEL - Parallel electron correlation program.
// Filename: libint.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifdef LIBINT_INTERFACE

#ifndef __SRC_RYSINT_LIBINT_H
#define __SRC_RYSINT_LIBINT_H

#include <src/integral/rys/rysintegral.h>

namespace bagel {

class Libint : public RysIntegral<double> {
  protected:
    void root_weight(int) {}
    void compute_ssss(double) {}

  public:
    Libint(const std::array<std::shared_ptr<const Shell>,4>&, const double dum = 0.0, std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>());

    void compute() {}

    constexpr static int Nblocks() { return 1; }

    virtual void allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) { assert(false); }
};

}

#endif

#endif

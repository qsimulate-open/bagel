//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI.h
// Copyright (C) 2015 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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


#ifndef __SRC_SMITH_RELMRCI_H
#define __SRC_SMITH_RELMRCI_H

#include <iostream>
#include <tuple>
#include <iomanip>
#include <src/smith/spinfreebase.h>
#include <src/smith/futuretensor.h>
#include <src/scf/hf/fock.h>
#include <src/util/f77.h>
#include <src/smith/queue.h>
#include <src/smith/multitensor.h>
#include <src/smith/smith_info.h>

namespace bagel {
namespace SMITH {
namespace RelMRCI{

class RelMRCI : public SpinFreeMethod<std::complex<double>> {
  protected:

  public:
    RelMRCI(std::shared_ptr<const SMITH_Info<std::complex<double>>> ref) : SpinFreeMethod(ref) { }
    ~RelMRCI() {}

    void solve() { }
    void solve_deriv() { }

    double accumulate(std::shared_ptr<Queue> queue) {
      return 0.0;
    }

};

}
}
}
#endif


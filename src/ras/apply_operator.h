//
// BAGEL - Parallel electron correlation program.
// Filename: ras/apply_operator.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __SRC_RAS_APPLY_OPERATOR_H
#define __SRC_RAS_APPLY_OPERATOR_H

#include <src/ciutil/citraits.h>
#include <src/ras/civector.h>
#include <src/asd/gamma_sq.h>

namespace bagel {

/// helper class used to apply second quantization operators used in plays like ASD_DMRG and similar in spirit to FormSigmaRAS
class ApplyOperator {
  private:
    /// takes a potential target string, returns whether it is valid and flips the proper bit in the process
    bool test_and_set(std::bitset<nbit__>& bit, const bool action, const int orbital) {
      bool out = action ? bit[orbital] : !bit[orbital];
      action ? bit.reset(orbital) : bit.set(orbital);
      return out;
    }

  public:
    ApplyOperator() {}

    void operator()(const RASCivecView source, RASCivecView target, std::vector<GammaSQ> operations, std::vector<int> orbitals) const;
};

}

#endif

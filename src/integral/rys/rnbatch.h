//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rnbatch.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_RYSINT_RNBATCH_H
#define __SRC_RYSINT_RNBATCH_H

#include <src/integral/rys/coulombbatch_energy.h>

namespace bagel {

class RnBatch: public CoulombBatch_energy {
  protected:
    int get_max_rterms() const {
      int out = 0;
      for (auto& aiter : mol_->atoms())
        if (aiter->use_ecp_basis())
          if (aiter->ecp_parameters()->shell_maxl_ecp()->ecp_exponents().size() > out)
            out = aiter->ecp_parameters()->shell_maxl_ecp()->ecp_exponents().size();
      return out;
    }

  public:
    RnBatch(const std::array<std::shared_ptr<const Shell>,2>& _info,
            const std::shared_ptr<const Molecule> mol, std::shared_ptr<StackMem> stack = nullptr)
      : CoulombBatch_energy (_info, mol, stack) {

       max_rterms_ = get_max_rterms();
       this->allocate_arrays(primsize_ * natom_ * max_rterms_);
    }

    ~RnBatch() { }

};

}

#endif

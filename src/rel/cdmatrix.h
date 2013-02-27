//
// BAGEL - Parallel electron correlation program.
// Filename: cdmatrix.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __src_rel_cdmatrix_h
#define __src_rel_cdmatrix_h

#include <cassert>
#include <src/util/matrix.h>
#include <src/util/zmatrix.h>
#include <src/rel/dfhalfcomplex.h>
#include <src/rel/reldfbase.h>

namespace bagel {

class CDMatrix : public ZMatrix {
  protected:
    const int comp_;

  public:
    CDMatrix(std::shared_ptr<DFHalfComplex> dfhc, std::shared_ptr<ABcases> abc, std::array<std::shared_ptr<const Matrix>, 4> trcoeff,
             std::array<std::shared_ptr<const Matrix>, 4> ticoeff, std::shared_ptr<const Matrix> dat2) : ZMatrix(
      *dfhc->get_real()->compute_cd(trcoeff[abc->basis(1)], dat2, true)+*dfhc->get_imag()->compute_cd(ticoeff[abc->basis(1)], dat2, true),
      *dfhc->get_real()->compute_cd(ticoeff[abc->basis(1)], dat2, true)-*dfhc->get_imag()->compute_cd(trcoeff[abc->basis(1)], dat2, true)), comp_(abc->comp()) {

      *this *= abc->fac();
    }

    CDMatrix(const ZMatrix& o, const int comp) : ZMatrix(o), comp_(comp) { }

    // multiply cd and breit2index for use in Jop in dfock.cc
    std::list<std::shared_ptr<const CDMatrix>> compute_breit_cd(std::list<std::shared_ptr<Breit2Index>>& b) const {
      std::list<std::shared_ptr<const CDMatrix>> out;
      for (auto i : b) {
        if (i->index().second == comp_)
          out.push_back(std::shared_ptr<CDMatrix>(new CDMatrix(*i->j_term() * *this, i->index().first)));
      }
      return out;
    }

    const int comp() const { return comp_; }
};

}

#endif

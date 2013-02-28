//
// BAGEL - Parallel electron correlation program.
// Filename: cdmatrix_drv.h
// Copyright (C) 2013 Matthew Kelley
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#ifndef __SRC_REL_CDMATRIX_DRV_H
#define __SRC_REL_CDMATRIX_DRV_H

#include <cassert>
#include <src/rel/dfhalfcomplex.h>
#include <src/rel/cdmatrix.h>

namespace bagel {

class CDMatrix_drv : public CDMatrix {
  protected:

  public:
    CDMatrix_drv(std::shared_ptr<DFHalfComplex> dfhc, std::shared_ptr<ABcases> abc, std::array<std::shared_ptr<const Matrix>, 4> trcoeff,
                 std::array<std::shared_ptr<const Matrix>, 4> ticoeff, std::shared_ptr<const Matrix> dat2)
    : CDMatrix(ZMatrix(*dfhc->get_real()->compute_cd(trcoeff[abc->basis(1)], dat2, true)+*dfhc->get_imag()->compute_cd(ticoeff[abc->basis(1)], dat2, true),
                       *dfhc->get_real()->compute_cd(ticoeff[abc->basis(1)], dat2, true)-*dfhc->get_imag()->compute_cd(trcoeff[abc->basis(1)], dat2, true)),
               abc->comp()) {

      *this *= abc->fac();
    }

};

}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: cdmatrix_london.cc
// Copyright (C) 2014 Toru Shiozaki
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


#include <src/london/cdmatrix_london.h>

using namespace std;
using namespace bagel;


// Also, C and D matrices are either real (for Coulomb) or purely imaginary (for Gaunt and Breit) due to symmetry. We are not taking advantage of it.

CDMatrix_London::CDMatrix_London(shared_ptr<const RelDFHalf_London> dfhc, shared_ptr<const SpinorInfo> abc, array<shared_ptr<const ZMatrix>, 4> trcoeff,
                                 array<shared_ptr<const ZMatrix>, 4> ticoeff, shared_ptr<const ZMatrix> dat2, const bool onlyonce)
   // TODO This is an absurd, ugly patch...  and will be incorrect if the two parts are not real.
 : ZMatrix(*(*dfhc->get_real()->compute_cd(trcoeff[abc->basis(1)], dat2, onlyonce)-*dfhc->get_imag()->compute_cd(ticoeff[abc->basis(1)], dat2, onlyonce)).get_real_part(),
           *(*dfhc->get_real()->compute_cd(ticoeff[abc->basis(1)], dat2, onlyonce)+*dfhc->get_imag()->compute_cd(trcoeff[abc->basis(1)], dat2, onlyonce)).get_real_part()),
   alpha_comp_(abc->alpha_comp()) {

  *this *= abc->fac(dfhc->cartesian());
}

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relcdmatrix.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/df/relcdmatrix.h>

using namespace std;
using namespace bagel;

// Also, C and D matrices are either real (for Coulomb) or purely imaginary (for Gaunt and Breit) due to symmetry. We are not taking advantage of it.

RelCDMatrix::RelCDMatrix(shared_ptr<const RelDFHalf> dfhc, shared_ptr<const SpinorInfo> abc, array<shared_ptr<const Matrix>, 4> trcoeff,
                         array<shared_ptr<const Matrix>, 4> ticoeff, shared_ptr<const Matrix> dat2, const int number_of_j)
 : ZVectorB(*dfhc->get_real()->compute_cd(trcoeff[abc->basis(1)], dat2, number_of_j)-*dfhc->get_imag()->compute_cd(ticoeff[abc->basis(1)], dat2, number_of_j),
            *dfhc->get_real()->compute_cd(ticoeff[abc->basis(1)], dat2, number_of_j)+*dfhc->get_imag()->compute_cd(trcoeff[abc->basis(1)], dat2, number_of_j)),
   alpha_comp_(abc->alpha_comp()) {

  btas::scal(abc->fac(dfhc->cartesian()), *this);

}


// BAGEL - Parallel electron correlation program.
// Filename: dimer_levelshift.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#ifndef __BAGEL_DIMER_DIMER_LEVELSHIFT_H
#define __BAGEL_DIMER_DIMER_LEVELSHIFT_H

#include <src/scf/coeff.h>
#include <src/scf/levelshift.h>
#include <src/dimer/dimer.h>

/************************************************************************************
* Class of functions to level shift the Fock operator in SCF. Matrices are shifted  *
*   using the shift(coeff) function, where coeff should be expressed in the         *
*   molecular basis                                                                 *
************************************************************************************/

namespace bagel {

// Dimer level shift
class ShiftDimer : public LevelShift<Matrix> {
  protected:
    std::shared_ptr<const Matrix> subspace_;
    std::shared_ptr<const Matrix> subspace_projector_;
    std::shared_ptr<const Matrix> S_;

    std::pair<int,int> nbasis_;
    int dimerbasis_;

  public:
    // monomer_coeff should be monomer coefficients projected onto dimer orbitals
    ShiftDimer(std::shared_ptr<const Dimer> dimer, const double shift_parameter);

    void shift(Matrix& fock_mo, std::shared_ptr<const Coeff> coeff);

    void print_mo_data(std::shared_ptr<const Coeff> coeff);
};

}

#endif


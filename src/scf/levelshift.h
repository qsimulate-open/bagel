
// BAGEL - Parallel electron correlation program.
// Filename: levelshift.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __BAGEL_SCF_LEVELSHIFT_H
#define __BAGEL_SCF_LEVELSHIFT_H

#include <array>
#include <vector>
#include <string>
#include <memory>

#include <src/scf/coeff.h>

/************************************************************************************
* Class of functions to level shift the Fock operator in SCF. Matrices are shifted  *
*   using the shift(coeff) function, where coeff should be expressed in the         *
*   molecular basis                                                                 *
************************************************************************************/

namespace bagel {

// This is the base class. It literally does nothing.
class LevelShift {
  protected:
    int nocc_;
    double shift_parameter_;

  public:
    LevelShift() {};
    LevelShift(const int nocc, const double shift_parameter) : nocc_(nocc), shift_parameter_(shift_parameter) {};

    virtual void shift(Matrix1e& coeff) {};
};

// Standard shifting of virtual orbitals
class ShiftVirtual : public LevelShift {
  protected:

  public:
    ShiftVirtual(const int nocc, const int shift_parameter);

    void shift(Matrix1e& coeff) override;
};

}

#endif

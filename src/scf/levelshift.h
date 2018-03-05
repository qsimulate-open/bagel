//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: levelshift.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __BAGEL_SCF_LEVELSHIFT_H
#define __BAGEL_SCF_LEVELSHIFT_H

#include <src/wfn/coeff.h>

/************************************************************************************
* Class of functions to level shift the Fock operator in SCF. Matrices are shifted  *
*   using the shift(coeff) function, where coeff should be expressed in the         *
*   molecular basis                                                                 *
************************************************************************************/

namespace bagel {

// This is the base class. It literally does nothing.
template<class Type>
class LevelShift {
  protected:
    int nocc_;
    double shift_parameter_;

  public:
    LevelShift() {};
    LevelShift(const int nocc, const double shift_parameter) : nocc_(nocc), shift_parameter_(shift_parameter) {};

    virtual void shift(Type& coeff) { }
};

// Standard shifting of virtual orbitals
template<class Type>
class ShiftVirtual : public LevelShift<Type> {
  protected:

  public:
    ShiftVirtual(const int nocc, const double shift_parameter) : LevelShift<Type>(nocc, shift_parameter) { }

    void shift(Type& coeff) override {
      coeff.add_diag(this->shift_parameter_, this->nocc_, coeff.ndim());
    }
};

}

#endif

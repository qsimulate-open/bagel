//
// BAGEL - Parallel electron correlation program.
// Filename: constant.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU Theory Group
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


#ifndef __SRC_UTIL_CONSTANTS_H
#define __SRC_UTIL_CONSTANTS_H

#include <cmath>
#include <chrono>
#include <stddef.h>

namespace bagel {

/************************************************************
*  Internal constants                                       *
************************************************************/
static constexpr int NDEBUG_PRINT = 1;
static constexpr int ANG_HRR_END = 7;
static constexpr int ANG_VRR_END = 14;
static constexpr int RYS_MAX = 13;
static constexpr double PRIM_SCREEN_THRESH = 1.0e-12;

/************************************************************
*  Fundamental Physical/Mathematical constants              *
************************************************************/
static const double pi__ = std::atan(1.0)*4.0;
static const double rad2deg__ = 180.0 / pi__;
static constexpr double c__ = 137.035999074;               // CODATA 2010 inverse fine-structure constant
static constexpr double au2kilogram__ = 9.10938291e-31;    // CODATA 2010 electron rest mass
static constexpr double au2coulomb__ = 1.602176565e-19;    // CODATA 2010 elementary charge
static constexpr double au2meter__ = 5.2917721092e-11;     // CODATA 2010 Bohr radius
static constexpr double avogadro__ = 6.02214129e23;        // CODATA 2010 Avogadro constant

/************************************************************
*  Derived unit conversions                                 *
************************************************************/
static const double au2second__ = c__ * au2meter__ / 2.99792458e8;
static const double au2angstrom__ = au2meter__ * 1.0e10;
static const double au2joule__ = au2kilogram__ * std::pow(au2meter__ / au2second__, 2);
static const double au2kjmol__ = au2joule__ * avogadro__ / 1.0e3;
static const double au2eV__ = au2kilogram__ * au2meter__ * au2meter__ / au2second__ / au2second__ / au2coulomb__;
static const double au2tesla__ = au2kilogram__ / au2coulomb__ / au2second__;

/************************************************************
*  Numerical constants                                      *
************************************************************/
static constexpr double numerical_zero__ = 1.0e-15;
static constexpr unsigned int nbit__ = 64;

/************************************************************
*  Misc constants                                           *
************************************************************/
static constexpr double schwarz_thresh__ = 1.0e-12;  // TODO input

/************************************************************
*  MPI parameters                                           *
************************************************************/
static constexpr size_t probe_key__  = (1 << 20);
static constexpr size_t probe_key2__ = (1 << 26);
static constexpr size_t pool_size__ = 100;
static constexpr std::chrono::microseconds sleeptime__ = std::chrono::microseconds(100);

}

#endif


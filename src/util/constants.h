//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: constant.h
// Copyright (C) 2012 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU Theory Group
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
static constexpr int RYS_MAX = 21;
static constexpr double PRIM_SCREEN_THRESH = 1.0e-12;

static constexpr int nucleus_blocksize__ = 500;            // maximum number of point charges used in a single batch of nuclear attraction integrals

/************************************************************
*  Fundamental Physical/Mathematical constants              *
************************************************************/
static const double pi__ = std::atan(1.0)*4.0;
static const double rad2deg__ = 180.0 / pi__;
static constexpr double c__ = 137.035999139;               // CODATA 2014 inverse fine-structure constant
static constexpr double au2kilogram__ = 9.10938356e-31;    // CODATA 2014 electron rest mass
static constexpr double au2coulomb__ = 1.6021766208e-19;   // CODATA 2014 elementary charge
static constexpr double au2meter__ = 5.2917721067e-11;     // CODATA 2014 Bohr radius
static constexpr double avogadro__ = 6.022140857e23;       // CODATA 2014 Avogadro constant
static constexpr double g_elec__ = 2.00231930436182;       // Absolute value of CODATA 2014 electron g factor

/************************************************************
*  Derived unit conversions                                 *
************************************************************/
static const double au2second__ = c__ * au2meter__ / 2.99792458e8;
static const double au2angstrom__ = au2meter__ * 1.0e10;
static const double au2joule__ = au2kilogram__ * std::pow(au2meter__ / au2second__, 2);
static const double au2kjmol__ = au2joule__ * avogadro__ / 1.0e3;
static const double au2eV__ = au2kilogram__ * au2meter__ * au2meter__ / au2second__ / au2second__ / au2coulomb__;
static const double au2tesla__ = au2kilogram__ / au2coulomb__ / au2second__;
static const double au2wavenumber__ = 1.0 / (2.0 * pi__ * c__ * au2meter__ * 100.0);

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


//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: constant.h
// Copyright (C) 2012 Toru Shiozaki
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
static constexpr int ANG_HRR_END = 8;
static constexpr int ANG_VRR_END = 16;
static constexpr int RYS_MAX = 21;
static constexpr double PRIM_SCREEN_THRESH = 1.0e-12;

static constexpr int nucleus_blocksize__ = 500;            // maximum number of point charges used in a single batch of nuclear attraction integrals

/************************************************************
*  Fundamental Physical/Mathematical constants              *
************************************************************/
static const double pi__ = std::atan(1.0)*4.0;
static const double rad2deg__ = 180.0 / pi__;
static constexpr double c__ = 137.035999139;               // CODATA 2014 inverse fine-structure constant
static constexpr double csi__ = 2.99792458e8;              // CODATA 2014 speed of light in vacuum
static constexpr double au2kilogram__ = 9.10938356e-31;    // CODATA 2014 electron rest mass
static constexpr double au2coulomb__ = 1.6021766208e-19;   // CODATA 2014 elementary charge
static constexpr double au2meter__ = 5.2917721067e-11;     // CODATA 2014 Bohr radius
static constexpr double avogadro__ = 6.022140857e23;       // CODATA 2014 Avogadro constant
static constexpr double g_elec__ = 2.00231930436182;       // Absolute value of CODATA 2014 electron g factor
static constexpr double amu2kilogram__ = 1.660539040e-27;  // COTDATA 2014 atomic mass unit-kilogram relationship
static constexpr double kcal2kj__ = 4.184;                 // Absolute definition of thermochemical calorie

/************************************************************
*  Derived unit conversions                                 *
************************************************************/
static const double au2second__ = c__ * au2meter__ / csi__;
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

}

#endif


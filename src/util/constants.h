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


#ifndef __src_util_constants_h
#define __src_util_constants_h

#include <cmath>

namespace bagel {

/************************************************************
*  Physical/Mathematical constants                          *
************************************************************/
static const double ang2bohr__ = 1.889725989;
static const double pi__ = std::atan(1.0)*4.0;
static const double rad2deg__ = 180.0 / pi__;

/************************************************************
*  Numerical constants                                      *
************************************************************/
static const double numerical_zero__ = 1.0e-15;
static const unsigned int large__ = 32;
static const unsigned int nbit__ = 32;

}

#endif


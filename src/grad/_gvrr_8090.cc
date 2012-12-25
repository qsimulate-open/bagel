//
// BAGEL - Parallel electron correlation program.
// Filename: _gvrr_8090.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/rysint/_vrr.h>
#include <src/grad/gvrrlist.h>

using namespace bagel;

// returns double array of length 810
void GVRRList::_gvrr_8090(double* data_, const double* C00, const double* D00, const double* B00, const double* B01, const double* B10) {
  vrr<8,9,9>(data_, C00, D00, B00, B01, B10);
}


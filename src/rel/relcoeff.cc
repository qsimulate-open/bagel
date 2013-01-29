//
// BAGEL - Parallel electron correlation program.
// Filename: relcoeff.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/rel/relcoeff.h>

using namespace std;
using namespace bagel;

RelCoeff::RelCoeff(const ZMatrix& inp) : ZMatrix(inp) {
}

array<shared_ptr<const ZMatrix>, 4> RelCoeff::split(const int n, const int m, const int msize) {
  array<shared_ptr<const ZMatrix>, 4> out;
  for (int i = 0; i != 4; ++i) {
    out[i] = shared_ptr<const ZMatrix>(new ZMatrix(*this->get_submatrix(i*n, m, n, msize))); 
  }

  return out;
}



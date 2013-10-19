#if 0
//
// BAGEL - Parallel electron correlation program.
// Filename: comperibatch.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Ryan Reynolds <rreynoldschem@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/integral/comprys/comperibatch.h>
#include <complex>

using namespace std;
using namespace bagel;


CompERIBatch::CompERIBatch(const array<shared_ptr<const Shell>,4>& _info, const double max_density, const double dummy, const bool dum,
                   shared_ptr<StackMem> stack) :  ERIBatch_Base<complex<double>>(_info, max_density, 0, 0, stack) {

#ifdef LIBINT_INTERFACE
  assert(false);
#endif
}


void CompERIBatch::compute() {
  throw logic_error("CompERIBatch::compute not yet implemented");
}
#endif

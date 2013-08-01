//
// BAGEL - Parallel electron correlation program.
// Filename: reldipole.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/rel/reldipole.h>
#include <src/molecule/dipolematrix.h>
#include <src/integral/os/dipolebatch.h>
#include <src/rel/small1e.h>

using namespace std;
using namespace bagel;

void RelDipole::compute() {
  // first calculate AO integrals.
  shared_ptr<const Matrix1eArray<3>> large = make_shared<const DipoleMatrix>(geom_);
  shared_ptr<const Matrix1eArray<12>> small = make_shared<const Small1e<DipoleBatch>>(geom_);

}

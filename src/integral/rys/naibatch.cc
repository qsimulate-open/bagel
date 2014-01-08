//
// BAGEL - Parallel electron correlation program.
// Filename: naibatch.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/integral/rys/naibatch.h>

using namespace std;
using namespace bagel;


NAIBatch::NAIBatch(const array<shared_ptr<const Shell>,2>& _info, const shared_ptr<const Molecule> mol, shared_ptr<StackMem> stack)
  : CoulombBatch_energy(_info, mol, stack) {
  const double integral_thresh = PRIM_SCREEN_THRESH;
  compute_ssss(integral_thresh);
  root_weight(primsize_*natom_);
}


NAIBatch::NAIBatch(const array<shared_ptr<const Shell>,2>& _info, const shared_ptr<const Molecule> mol, const int L, const double A)
  : CoulombBatch_energy (_info, mol, L, A) {
  const double integral_thresh = PRIM_SCREEN_THRESH;
  compute_ssss(integral_thresh);
  root_weight(primsize_*natom_);
}


//
// BAGEL - Parallel electron correlation program.
// Filename: multipolebatch_base.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/multipolebatch_base.h>

using namespace std;
using namespace bagel;

MultipoleBatch_base::MultipoleBatch_base(const array<shared_ptr<const Shell>,2>& sh, const shared_ptr<const Atom> a,
                                         shared_ptr<StackMem> stack)
 : basisinfo_(sh), site_(a), spherical_(sh[0]->spherical()) {
  assert(spherical_ = sh[1]->spherical());

  if (stack == nullptr) {
    stack_ = resources__->get();
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }

}


void MultipoleBatch_base::compute_ss(const double thr) {


}

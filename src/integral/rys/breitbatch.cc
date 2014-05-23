//
// BAGEL - Parallel electron correlation program.
// Filename: breitbatch.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/integral/rys/breitbatch.h>
#include <src/integral/rys/breitrootlist.h>
#include <src/integral/rys/spin2rootlist.h>

using namespace std;
using namespace bagel;

void BreitBatch_base::root_weight(const int ps) {
  if (breit_ == 1) {
    breitroot__.root(rank_, T_, roots_, weights_, ps);
  } else if (breit_ == 2) {
    spin2root__.root(rank_, T_, roots_, weights_, ps);
  } else {
    assert(false);
  }
}


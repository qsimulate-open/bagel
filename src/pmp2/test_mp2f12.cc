//
// Newint - Parallel electron correlation program.
// Filename: test_mp2f12.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <src/pscf/pgeometry.h>
#include <src/pscf/poverlap.h>
#include <src/pscf/pscf.h>
#include <src/pscf/pscf_disk.h>
#include <src/pmp2/pmp2.h>

#include <src/stackmem.h>

#include <src/util/input.h>

using namespace std;
extern  StackMem* stack;

void test_mp2f12() {
  StackMem* a = new StackMem(static_cast<size_t>(100000000LU));
  stack = a;

  shared_ptr<PGeometry> pgeom(new PGeometry("oldinp", 0));
  shared_ptr<PSCF_DISK> pscf(new PSCF_DISK(pgeom));
  pscf->compute();

  PMP2 pmp2(pscf->geom(), pscf->coeff(), pscf->eig(), pscf->ao_eri(), false);
  pmp2.compute();


};


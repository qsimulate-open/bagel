//
// BAGEL - Parallel electron correlation program.
// Filename: pscf_disk.cc
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


#include <src/pscf/pscf_disk.h>
#include <src/integral/rys/eribatch.h>
#include <src/pscf/pcompfile.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstring>

using namespace std;
using namespace bagel;

typedef shared_ptr<const Atom> RefAtom;
typedef shared_ptr<Shell> RefShell;
typedef shared_ptr<PGeometry> RefPGeometry;
typedef shared_ptr<PFock> RefPFock;
typedef shared_ptr<PTildeX> RefPTildeX;
typedef shared_ptr<PCoeff> RefPCoeff;
typedef shared_ptr<PMatrix1e> RefPMatrix1e;

PSCF_DISK::PSCF_DISK(const RefPGeometry g) : PSCF(g) {

  store_ERI();

}


PSCF_DISK::PSCF_DISK(const RefPGeometry g, RefPMatrix1e a) : PSCF(g, a) {

  store_ERI();

}

PSCF_DISK::~PSCF_DISK() {

}


void PSCF_DISK::store_ERI() {

  shared_ptr<PCompFile<ERIBatch>> tmpae(new PCompFile<ERIBatch>(geom_, 0.0, false, "ERI OBS"));
  ao_eri_ = tmpae;
  ao_eri_->store_integrals();
  ao_eri_->reopen_with_inout();

}


//
// BAGEL - Parallel electron correlation program.
// Filename: erirootlist.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/rysint/erirootlist.h>
#include <src/rysint/f77.h>

using namespace bagel;

ERIRootList::ERIRootList() {

  rfunc[1] = &eriroot1_;
  rfunc[2] = &eriroot2_;
  rfunc[3] = &eriroot3_;
  rfunc[4] = &eriroot4_;
  rfunc[5] = &eriroot5_;
  rfunc[6] = &eriroot6_;
  rfunc[7] = &eriroot7_;
  rfunc[8] = &eriroot8_;
  rfunc[9] = &eriroot9_;
  rfunc[10] = &eriroot10_;
  rfunc[11] = &eriroot11_;
  rfunc[12] = &eriroot12_;
}


ERIRootList::~ERIRootList() {

}

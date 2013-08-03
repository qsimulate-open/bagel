//
// BAGEL - Parallel electron correlation program.
// Filename: smith.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/smith/mp2_ref.h>
#include <src/smith/storage.h>

#include <src/smith/MP2.h>
#include <src/smith/CAS_all_active.h>

using namespace bagel;
using namespace bagel::SMITH;
using namespace std;

// called from the main function
void smith_test(shared_ptr<Reference> r) {
  {
    MP2_Ref<Storage_Incore> mp2(r);
    mp2.solve();
  }

  {
    MP2::MP2<Storage_Incore> _mp2(r);
    _mp2.solve();
  }
  // in order to see if it compiles
  if (false) {
    CAS_all_active::CAS_all_active<Storage_Incore> cas(r);
  }
}



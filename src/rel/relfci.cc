//
// BAGEL - Parallel electron correlation program.
// Filename: relfci.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Matthew Kelley matthewkelley2017@u.northwestern.edu
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

#include <string>
#include <vector>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include <src/rel/relfci.h>
#include <src/fci/space.h>
#include <src/rysint/eribatch.h>
#include <src/util/combination.hpp>
#include <src/util/constants.h>
#include <src/util/davidson.h>

using namespace std;
using namespace bagel;

RelFCI::RelFCI(std::multimap<std::string, std::string> idat, shared_ptr<const RelReference> r)
 : idata_(idat), ref_(r), geom_(r->geom()) {
}

void RelFCI::compute() {
  cout << "Nothing here yet... " << endl;

}


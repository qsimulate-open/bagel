//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/construct_asd_oo.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
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

#include <src/asd/orbital/construct_asd_oo.h>
#include <src/asd/orbital/bfgs.h>

using namespace std;
using namespace bagel;

namespace bagel {

shared_ptr<ASD_OO> construct_ASD_OO(shared_ptr<const PTree> itree, shared_ptr<Dimer> dimer) {
  string algorithm = itree->get<string>("algorithm", "");
  string method = itree->get_child_optional("asd")->get<string>("method");

  shared_ptr<ASD_OO> out;

  if (method == "cas" && algorithm == "bfgs") {
   out = make_shared<ASD_BFGS>(itree, dimer);
  } else if (method == "ras" && algorithm == "bfgs") {
   out = make_shared<ASD_BFGS>(itree, dimer);
  } else {
    throw runtime_error("unsupported combination of ASD orbital optimization algorithm specified: " + method + "+" + algorithm);
  }

  return out;
}

}


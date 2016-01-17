//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd/orbital/construct_asd_orbopt.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/asd/orbital/construct_asd_orbopt.h>
#include <src/asd/orbital/asd_bfgs.h>

using namespace std;
using namespace bagel;

namespace bagel {

shared_ptr<ASD_OrbOpt> construct_ASD_OrbOpt(shared_ptr<const PTree> itree, shared_ptr<Dimer> dimer) {
  string algorithm = itree->get<string>("algorithm", "");
  string method = itree->get_child_optional("asd")->get<string>("method");

  shared_ptr<ASD_OrbOpt> out;

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


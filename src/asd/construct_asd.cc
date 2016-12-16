//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: construct_asd.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/asd_cas.h>
#include <src/asd/asd_ras.h>

using namespace std;
using namespace bagel;

namespace bagel {

shared_ptr<ASD_base> construct_ASD(shared_ptr<const PTree> itree, shared_ptr<Dimer> dimer, bool rdm) {
  shared_ptr<ASD_base> out;
  string method = itree->get<string>("method");
  string variant = itree->get<string>("variant", "local");

  if (method == "cas" || method == "fci") {
    if (variant == "local") {
      shared_ptr<DimerCAS> cispace = dimer->compute_cispace<CASDvec>(itree);
      out = make_shared<ASD_CAS>(itree, dimer, cispace, rdm);
    } else {
      throw logic_error("Unrecognized variant of CAS ASD");
    }
  } else if (method == "ras") {
    if (variant == "local") {
      shared_ptr<DimerRAS> cispace = dimer->compute_rcispace<RASDvec>(itree);
      out = make_shared<ASD_RAS>(itree, dimer, cispace, rdm);
    } else {
      throw logic_error("Unrecognized variant of RAS ASD");
    }
  } else {
    throw runtime_error("Unrecognized method for ASD");
  }

  return out;
}

}

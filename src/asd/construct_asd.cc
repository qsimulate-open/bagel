//
// BAGEL - Parallel electron correlation program.
// Filename: construct_asd.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/asd_cas.h>
#include <src/asd/asd_distcas.h>
#include <src/asd/asd_ras.h>
#include <src/asd/asd_distras.h>

using namespace std;
using namespace bagel;

namespace bagel {

shared_ptr<ASD_base> construct_ASD(shared_ptr<const PTree> itree, shared_ptr<Dimer> dimer) {
  shared_ptr<ASD_base> out;
  string method = itree->get<string>("method");
  string variant = itree->get<string>("variant", "local");

  if (method == "cas" || method == "fci") {
    if (variant == "local") {
      shared_ptr<DimerCAS> cispace = dimer->compute_cispace<Dvec>(itree);
      out = make_shared<ASD_CAS>(itree, dimer, cispace);
    }
    else if (variant == "dist" || variant == "parallel") {
      shared_ptr<DimerDistCAS> cispace = dimer->compute_cispace<DistDvec>(itree);
      out = make_shared<ASD_DistCAS>(itree, dimer, cispace);
    }
    else {
      throw logic_error("Unrecognized variant of CAS ASD");
    }
  }
  else if (method == "ras") {
    if (variant == "local") {
      shared_ptr<DimerRAS> cispace = dimer->compute_rcispace<RASDvec>(itree);
      out = make_shared<ASD_RAS>(itree, dimer, cispace);
    }
    else if (variant == "dist" || variant == "parallel") {
      shared_ptr<DimerDistRAS> cispace = dimer->compute_rcispace<DistRASDvec>(itree);
      out = make_shared<ASD_DistRAS>(itree, dimer, cispace);
    }
    else {
      throw logic_error("Unrecognized variant of RAS ASD");
    }
  }
  else {
    throw runtime_error("Unrecognized method for ASD");
  }

  return out;
}

}

//
// BAGEL - Parallel electron correlation program.
// Filename: construct_meh.cc
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

#include <src/meh/meh_cas.h>
#include <src/meh/meh_distcas.h>
#include <src/meh/meh_ras.h>
#include <src/meh/meh_distras.h>

using namespace std;
using namespace bagel;

namespace bagel {

shared_ptr<MEH_base> construct_MEH(shared_ptr<const PTree> itree, shared_ptr<Dimer> dimer) {
  shared_ptr<MEH_base> out;
  string method = itree->get<string>("method");
  string variant = itree->get<string>("variant", "local");

  if (method == "cas" || method == "fci") {
    if (variant == "local") {
      shared_ptr<DimerCAS> cispace = dimer->compute_cispace<Dvec>(itree);
      out = make_shared<MEH_CAS>(itree, dimer, cispace);
    }
    else if (variant == "dist" || variant == "parallel") {
      shared_ptr<DimerDistCAS> cispace = dimer->compute_cispace<DistDvec>(itree);
      out = make_shared<MEH_DistCAS>(itree, dimer, cispace);
    }
    else {
      throw logic_error("Unrecognized variant of CAS MEH");
    }
  }
  else if (method == "ras") {
    if (variant == "local") {
      shared_ptr<DimerRAS> cispace = dimer->compute_rcispace<RASDvec>(itree);
      out = make_shared<MEH_RAS>(itree, dimer, cispace);
    }
    else if (variant == "dist" || variant == "parallel") {
      shared_ptr<DimerDistRAS> cispace = dimer->compute_rcispace<DistRASDvec>(itree);
      out = make_shared<MEH_DistRAS>(itree, dimer, cispace);
    }
    else {
      throw logic_error("Unrecognized variant of RAS MEH");
    }
  }
  else {
    throw runtime_error("Unrecognized method for MEH");
  }

  return out;
}

}

//
// BAGEL - Parallel electron correlation program.
// Filename: construct_method.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/scf/rohf.h>
#include <src/ks/ks.h>
#include <src/scf/soscf.h>
#include <src/fci/distfci.h>
#include <src/fci/harrison.h>
#include <src/fci/knowles.h>
#include <src/ras/rasci.h>
#include <src/ras/distrasci.h>
#include <src/zfci/zharrison.h>
#include <src/casscf/superci.h>
#include <src/casscf/werner.h>
#include <src/casscf/casbfgs.h>
#include <src/zcasscf/zcasscf.h>
#include <src/rel/dirac.h>
#include <src/rel/dmp2.h>
#include <src/mp2/mp2.h>
#include <src/smith/smith.h>
#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

namespace bagel {

shared_ptr<Method> construct_method(string title, shared_ptr<const PTree> itree, shared_ptr<const Geometry> geom,
                                                  shared_ptr<const Reference> ref) {
  shared_ptr<Method> out;
  if (title == "hf")          out = make_shared<SCF>(itree, geom, ref);
  else if (title == "ks")     out = make_shared<KS>(itree, geom, ref);
  else if (title == "uhf")    out = make_shared<UHF>(itree, geom, ref);
  else if (title == "rohf")   out = make_shared<ROHF>(itree, geom, ref);
  else if (title == "soscf")  out = make_shared<SOSCF>(itree, geom, ref);
  else if (title == "mp2")    out = make_shared<MP2>(itree, geom, ref);
  else if (title == "dhf")    out = make_shared<Dirac>(itree, geom, ref);
  else if (title == "dmp2")   out = make_shared<DMP2>(itree, geom, ref);
  else if (title == "smith")  out = make_shared<Smith>(itree, geom, ref);
  else if (title == "zfci")   out = make_shared<ZHarrison>(itree, geom, ref);
  else if (title == "ras") {
    const string algorithm = itree->get<string>("algorithm", "");
    if ( algorithm == "local" || algorithm == "" ) {
      out = make_shared<RASCI>(itree, geom, ref);
    }
#ifdef HAVE_MPI_H
    else if ( algorithm == "dist" || algorithm == "parallel" ) {
      out = make_shared<DistRASCI>(itree, geom, ref);
    }
#endif
    else
      throw runtime_error("unknown RASCI algorithm specified. " + algorithm);
  }
  else if (title == "fci") {
    const string algorithm = itree->get<string>("algorithm", "");
    const bool dokh = (algorithm == "" || algorithm == "auto") && geom->nele() > geom->nbasis();
    if (dokh || algorithm == "kh" || algorithm == "knowles" || algorithm == "handy") {
      out = make_shared<KnowlesHandy>(itree, geom, ref);
    } else if (algorithm == "hz" || algorithm == "harrison" || algorithm == "zarrabian" || algorithm == "") {
      out = make_shared<HarrisonZarrabian>(itree, geom, ref);
#ifdef HAVE_MPI_H
    } else if (algorithm == "parallel" || algorithm == "dist") {
      out = make_shared<DistFCI>(itree, geom, ref);
#endif
    } else
      throw runtime_error("unknown FCI algorithm specified. " + algorithm);
  }
  else if (title == "casscf") {
    string algorithm = itree->get<string>("algorithm", "");
    if (algorithm == "superci" || algorithm == "")
      out = make_shared<SuperCI>(itree, geom, ref);
    else if (algorithm == "werner" || algorithm == "knowles")
      out = make_shared<WernerKnowles>(itree, geom, ref);
    else if (algorithm == "bfgs")
      out = make_shared<CASBFGS>(itree, geom, ref);
    else
      throw runtime_error("unknown CASSCF algorithm specified: " + algorithm);
  }
  else if (title == "zcasscf") out = make_shared<ZCASSCF>(itree, geom, ref);

  return out;
}

}

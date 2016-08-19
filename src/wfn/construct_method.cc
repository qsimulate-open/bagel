//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: construct_method.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/scf/hf/rohf.h>
#include <src/scf/ks/ks.h>
#include <src/scf/sohf/soscf.h>
#include <src/scf/dhf/dirac.h>
#include <src/scf/giaohf/rhf_london.h>
#include <src/ci/fci/distfci.h>
#include <src/ci/fci/harrison.h>
#include <src/ci/fci/knowles.h>
#include <src/ci/ras/rasci.h>
#include <src/ci/ras/distrasci.h>
#include <src/ci/zfci/zharrison.h>
#include <src/pt2/nevpt2/nevpt2.h>
#include <src/pt2/mp2/mp2.h>
#include <src/pt2/dmp2/dmp2.h>
#include <src/multi/casscf/superci.h>
#include <src/multi/casscf/cashybrid.h>
#include <src/multi/casscf/casbfgs.h>
#include <src/multi/casscf/casnoopt.h>
#include <src/multi/zcasscf/zcasscf.h>
#include <src/multi/zcasscf/zcasbfgs.h>
#include <src/multi/zcasscf/zcashybrid.h>
#include <src/multi/zcasscf/zsuperci.h>
#include <src/multi/zcasscf/zcasnoopt.h>
#include <src/smith/smith.h>
#include <src/smith/caspt2grad.h>
#include <src/prop/current.h>
#include <src/prop/moprint.h>
#include <src/wfn/construct_method.h>

using namespace std;
using namespace bagel;

namespace bagel {

shared_ptr<Method> construct_method(string title, shared_ptr<const PTree> itree, shared_ptr<const Geometry> geom,
                                                  shared_ptr<const Reference> ref) {

  shared_ptr<Method> out;
  if (!geom->magnetism()) {
    if (title == "hf")            out = make_shared<RHF>(itree, geom, ref);
    else if (title == "ks")       out = make_shared<KS>(itree, geom, ref);
    else if (title == "uhf")      out = make_shared<UHF>(itree, geom, ref);
    else if (title == "rohf")     out = make_shared<ROHF>(itree, geom, ref);
    else if (title == "soscf")    out = make_shared<SOSCF>(itree, geom, ref);
    else if (title == "mp2")      out = make_shared<MP2>(itree, geom, ref);
    else if (title == "dhf")      out = make_shared<Dirac>(itree, geom, ref);
    else if (title == "dmp2")     out = make_shared<DMP2>(itree, geom, ref);
    else if (title == "smith")    out = make_shared<Smith>(itree, geom, ref);
    else if (title == "relsmith") out = make_shared<RelSmith>(itree, geom, ref);
    else if (title == "zfci")     out = make_shared<ZHarrison>(itree, geom, ref);
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
      else if (algorithm == "hybrid")
        out = make_shared<CASHybrid>(itree, geom, ref);
      else if (algorithm == "bfgs")
        out = make_shared<CASBFGS>(itree, geom, ref);
      else if (algorithm == "noopt")
        out = make_shared<CASNoopt>(itree, geom, ref);
      else
        throw runtime_error("unknown CASSCF algorithm specified: " + algorithm);
    }
    else if (title == "caspt2grad") {
      // TODO to be called from optimizer
      out = make_shared<CASPT2Grad>(itree, geom, ref);
    }
    else if (title == "nevpt2")  out = make_shared<NEVPT2<double>>(itree, geom, ref);
    else if (title == "dnevpt2") out = make_shared<NEVPT2<complex<double>>>(itree, geom, ref);
    else if (title == "zcasscf") {
      string algorithm = itree->get<string>("algorithm", "");
      if (algorithm == "superci" || algorithm == "")
        out = make_shared<ZSuperCI>(itree, geom, ref);
      else if (algorithm == "hybrid")
        out = make_shared<ZCASHybrid>(itree, geom, ref);
      else if (algorithm == "bfgs")
        out = make_shared<ZCASBFGS>(itree, geom, ref);
      else if (algorithm == "noopt")
        out = make_shared<ZCASNoopt>(itree, geom, ref);
      else
        cout << " Optimization algorithm " << algorithm << " is not compatible with ZCASSCF " << endl;
    } else if (title == "current")  throw runtime_error("Charge currents are only available when using a GIAO basis set reference.");
    else if (title == "moprint") out = make_shared<MOPrint>(itree, geom, ref);

  // now the versions to use with magnetic fields
  } else {
    if (title == "hf")              out = make_shared<RHF_London>(itree, geom, ref);
    else if (title == "dhf")        out = make_shared<Dirac>(itree, geom, ref);
    else if (title == "current")    out = make_shared<Current>(itree, geom, ref);
    else if (title == "zfci")       out = make_shared<ZHarrison>(itree, geom, ref);
    else if (title == "zcasscf") {
      string algorithm = itree->get<string>("algorithm", "");
      if (algorithm == "superci" || algorithm == "")
        out = make_shared<ZSuperCI>(itree, geom, ref);
      else if (algorithm == "hybrid")
        out = make_shared<ZCASHybrid>(itree, geom, ref);
      else if (algorithm == "bfgs")
        out = make_shared<ZCASBFGS>(itree, geom, ref);
      else
        cout << " Optimization algorithm " << algorithm << " is not compatible with ZCASSCF " << endl;
    } else if (title == "moprint") { out = make_shared<MOPrint>(itree, geom, ref);
    } else if (title == "molecule") {
    } else
      throw runtime_error(to_upper(title) + " method has not been implemented with an applied magnetic field.");
  }
  return out;
}

}

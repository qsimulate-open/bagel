//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: get_energy.cc
// Copyright (C) 2017 Toru Shiozaki
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
#include <src/ci/zfci/relfci.h>
#include <src/ci/zfci/fci_london.h>
#include <src/response/cis.h>
#include <src/pt2/nevpt2/nevpt2.h>
#include <src/pt2/mp2/mp2.h>
#include <src/pt2/dmp2/dmp2.h>
#include <src/multi/casscf/cassecond.h>
#include <src/multi/casscf/casnoopt.h>
#include <src/multi/zcasscf/zcassecond.h>
#include <src/multi/zcasscf/zcasnoopt.h>
#include <src/smith/smith.h>
#include <src/smith/caspt2grad.h>
#include <src/prop/current.h>
#include <src/prop/moprint.h>
#include <src/wfn/get_energy.h>

using namespace std;
using namespace bagel;

namespace bagel {

tuple<double, shared_ptr<const Reference>>
get_energy(const string title, shared_ptr<const PTree> itree, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref, const int target) {
  double out = 0.0;

#ifndef DISABLE_SERIALIZATION
  if (title == "continue") {
    shared_ptr<Method> m;
    IArchive archive(itree->get<string>("archive"));
    Method* ptr;
    archive >> ptr;
    m = shared_ptr<Method>(ptr);
    m->compute();
    ref = m->conv_to_ref();

    return tie(out, ref);
  }
#endif

  if (!geom || !geom->magnetism()) {
    if (title == "hf")           { auto m = make_shared<RHF>(itree, geom, ref);       m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
    else if (title == "ks")      { auto m = make_shared<KS>(itree, geom, ref);        m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
    else if (title == "uhf")     { auto m = make_shared<UHF>(itree, geom, ref);       m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
    else if (title == "rohf")    { auto m = make_shared<ROHF>(itree, geom, ref);      m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
    else if (title == "soscf")   { auto m = make_shared<SOSCF>(itree, geom, ref);     m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
    else if (title == "mp2")     { auto m = make_shared<MP2>(itree, geom, ref);       m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
    else if (title == "dhf")     { auto m = make_shared<Dirac>(itree, geom, ref);     m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
    else if (title == "dmp2")    { auto m = make_shared<DMP2>(itree, geom, ref);      m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
    else if (title == "cis")     { auto m = make_shared<CIS>(itree, geom, ref);       m->compute();   out = m->energy();                 ref = m->conv_to_ref(); }
#ifdef COMPILE_SMITH
    else if (title == "smith")   { auto m = make_shared<Smith>(itree, geom, ref);     m->compute();   out = m->algo()->energy(target);   ref = m->conv_to_ref(); }
    else if (title == "relsmith"){ auto m = make_shared<RelSmith>(itree, geom, ref);  m->compute();   out = m->algo()->energy(target);   ref = m->conv_to_ref(); }
#else
    else if (title == "smith" || title == "relsmith")   { throw runtime_error("SMITH module was not activated during compilation."); }
#endif
    else if (title == "zfci")    { auto m = make_shared<RelFCI>(itree, geom, ref);    m->compute();   out = m->energy(target);           ref = m->conv_to_ref(); }
    else if (title == "ras") {
      const string algorithm = itree->get<string>("algorithm", "");
      if ( algorithm == "local" || algorithm == "" ) { auto m = make_shared<RASCI>(itree, geom, ref); m->compute(); out = m->energy(target); ref = m->conv_to_ref(); }
#ifdef HAVE_MPI_H
      else if ( algorithm == "dist" || algorithm == "parallel" )
        throw logic_error("Parallel RASCI not implemented");
#endif
      else
        throw runtime_error("unknown RASCI algorithm specified. " + algorithm);
    }
    else if (title == "fci") {
      const string algorithm = itree->get<string>("algorithm", "");
      const bool dokh = (algorithm == "" || algorithm == "auto") && geom->nele() > geom->nbasis();
      if (dokh || algorithm == "kh" || algorithm == "knowles" || algorithm == "handy") {
        auto m = make_shared<KnowlesHandy>(itree, geom, ref);        m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
      } else if (algorithm == "hz" || algorithm == "harrison" || algorithm == "zarrabian" || algorithm == "") {
        auto m = make_shared<HarrisonZarrabian>(itree, geom, ref);   m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
#ifdef HAVE_MPI_H
      } else if (algorithm == "parallel" || algorithm == "dist") {
        auto m = make_shared<DistFCI>(itree, geom, ref);             m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
#endif
      } else
        throw runtime_error("unknown FCI algorithm specified. " + algorithm);
    }
    else if (title == "casscf") {
      string algorithm = itree->get<string>("algorithm", "");
      if (algorithm == "second" || algorithm == "") {
        auto m = make_shared<CASSecond>(itree, geom, ref);           m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
      } else if (algorithm == "noopt") {
        auto m = make_shared<CASNoopt>(itree, geom, ref);            m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
      } else
        throw runtime_error("unknown CASSCF algorithm specified: " + algorithm);
    }
    else if (title == "nevpt2")  { auto m = make_shared<NEVPT2<double>>(itree, geom, ref);          m->compute();   out = m->energy();    ref = m->conv_to_ref(); }
    else if (title == "dnevpt2") { auto m = make_shared<NEVPT2<complex<double>>>(itree, geom, ref); m->compute();   out = m->energy();    ref = m->conv_to_ref(); }
    else if (title == "zcasscf") {
      string algorithm = itree->get<string>("algorithm", "");
      if (algorithm == "second" || algorithm == "") {
        auto m = make_shared<ZCASSecond>(itree, geom, ref);          m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
      } else if (algorithm == "noopt") {
        auto m = make_shared<ZCASNoopt>(itree, geom, ref);           m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
      } else
        cout << " Optimization algorithm " << algorithm << " is not compatible with ZCASSCF " << endl;
    } else if (title == "current")  throw runtime_error("Charge currents are only available when using a GIAO basis set reference.");
    else if (title == "moprint") { auto m = make_shared<MOPrint>(itree, geom, ref);          m->compute();                         ref = m->conv_to_ref();
    } else if (title == "molecule") {
    } else throw runtime_error(to_upper(title) + " : unknown method");

  // now the versions to use with magnetic fields
  } else {
    if (title == "hf")           { auto m = make_shared<RHF_London>(itree, geom, ref);       m->compute();   out = m->energy();        ref = m->conv_to_ref(); }
    else if (title == "dhf")     { auto m = make_shared<Dirac>(itree, geom, ref);            m->compute();   out = m->energy();        ref = m->conv_to_ref(); }
    else if (title == "current") { auto m = make_shared<Current>(itree, geom, ref);          m->compute();                             ref = m->conv_to_ref(); }
    else if (title == "fci")     { auto m = make_shared<FCI_London>(itree, geom, ref);       m->compute();   out = m->energy(target);  ref = m->conv_to_ref(); }
    else if (title == "zfci")    { auto m = make_shared<RelFCI>(itree, geom, ref);           m->compute();   out = m->energy(target);  ref = m->conv_to_ref(); }
    else if (title == "zcasscf") {
      string algorithm = itree->get<string>("algorithm", "");
      if (algorithm == "second" || algorithm == "") {
        auto m = make_shared<ZCASSecond_London>(itree, geom, ref);   m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
      } else if (algorithm == "noopt") {
        auto m = make_shared<ZCASNoopt_London>(itree, geom, ref);    m->compute();   out = m->energy(target);   ref = m->conv_to_ref();
      } else
        cout << " Optimization algorithm " << algorithm << " is not compatible with ZCASSCF " << endl;
    } else if (title == "moprint") { auto m = make_shared<MOPrint>(itree, geom, ref);       m->compute();       ref = m->conv_to_ref();
    } else if (title == "molecule") {
    } else throw runtime_error(to_upper(title) + " method has not been implemented with an applied magnetic field.");
  }

  return tie(out, ref);
}

}

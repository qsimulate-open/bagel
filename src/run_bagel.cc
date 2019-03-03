//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: run_bagel.cc
// Copyright (C) 2009 Toru Shiozaki
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

#include <sstream>
#include <src/bagel.h>
#include <src/global.h>
#include <src/pt2/mp2/mp2grad.h>
#include <src/grad/force.h>
#include <src/grad/hess.h>
#include <src/opt/optimize.h>
#include <src/wfn/localization.h>
#include <src/asd/construct_asd.h>
#include <src/asd/orbital/construct_asd_orbopt.h>
#include <src/asd/dmrg/rasd.h>
#include <src/wfn/molden_to_ref.h>
#include <src/asd/multisite/multisite.h>
#include <src/util/archive.h>
#include <src/util/io/moldenout.h>

using namespace std;
using namespace bagel;

namespace bagel { namespace impl {
  void run_bagel_(shared_ptr<const PTree> idata);
} }


void bagel::run_bagel_from_input(const string& input) {
  static_variables();
  auto idata = make_shared<const PTree>(input);
  impl::run_bagel_(idata);
}


void bagel::run_bagel_from_json(const string& input) {
  static_variables();
  stringstream ss; ss << input;
  auto idata = make_shared<const PTree>(ss);
  impl::run_bagel_(idata);
}


void bagel::impl::run_bagel_(shared_ptr<const PTree> idata) {

  print_header();

  shared_ptr<const Geometry> geom;
  shared_ptr<const Reference> ref;
  shared_ptr<Dimer> dimer;

  map<string, shared_ptr<const void>> saved;
  bool dodf = true;
  bool dofmm = false;

  // timer for each method
  Timer timer(-1);

  shared_ptr<const PTree> keys = idata->get_child("bagel");

  for (auto& itree : *keys) {

    const string title = to_lower(itree->get<string>("title", ""));
    if (title.empty()) throw runtime_error("title is missing in one of the input blocks");

    if (title == "molecule") {
      geom = geom ? make_shared<Geometry>(*geom, itree) : make_shared<Geometry>(itree);
      if (itree->get<bool>("restart", false))
        ref.reset();
      if (ref) ref = ref->project_coeff(geom);
      if (!itree->get<string>("molden_file", "").empty())
        ref = molden_to_ref(geom, itree);
    } else {
      if (!geom) {
        if (title != "continue" && !(title == "relsmith" && itree->get<string>("method", "") == "continue") && title != "load_ref")
          throw runtime_error("molecule block is missing");
      } else {
        if (!itree->get<bool>("df",true)) dodf = false;
        if (dodf && !geom->df() && !geom->do_periodic_df() && !dofmm) throw runtime_error("It seems that DF basis was not specified in molecule block");
      }
    }

    if ((title == "smith" || title == "relsmith" || title == "fci") && ref == nullptr)
      if (itree->get<string>("method", "") != "continue")
        throw runtime_error(title + " needs a reference");

    if (title == "optimize") {

      auto opt = make_shared<Optimize>(itree, geom, ref);
      opt->compute();
      ref = opt->conv_to_ref();
      geom = opt->geometry();

    } else if (title == "forces" || title == "force" || title == "nacme" || title == "dgrad") {
      // "forces" does many force calculations at once (SA-CASSCF and MS-CASPT2)

      auto force = make_shared<Force>(itree, geom, ref);
      force->compute();
      ref = force->conv_to_ref();

    } else if (title == "hessian") {

      auto hess = make_shared<Hess>(itree, geom, ref);
      hess->compute();
      ref = hess->conv_to_ref();

#ifndef DISABLE_SERIALIZATION
    } else if (title == "load_ref") {
      const string name = itree->get<string>("file", "reference");
      if (name == "") throw runtime_error("Please provide a filename for the Reference object to be read.");
      IArchive archive(name);
      shared_ptr<Reference> ptr;
      archive >> ptr;
      ref = shared_ptr<Reference>(ptr);
      if (itree->get<bool>("continue_geom", true)) {
        cout << endl << "  Using the geometry in the archive file " << endl << endl;
        geom = ref->geom();
      } else {
        cout << endl << "  Using the coefficient projected to the input geometry " << endl << endl;
        ref = ref->project_coeff(geom);
      }
      if (itree->get<bool>("extract_average_rdms", false)) {
        vector<int> rdm_states = itree->get_vector<int>("rdm_state");
        ref = ref->extract_average_rdm(rdm_states);
      }

    } else if (title == "save_ref") {
      const string name = itree->get<string>("file", "reference");
      OArchive archive(name);
      archive << ref;
#endif
    } else if (title == "dimerize") { // dimerize forms the dimer object, does a scf calculation, and then localizes
      const string form = itree->get<string>("form", "displace");
      if (form == "d" || form == "disp" || form == "displace") {
        if (static_cast<bool>(ref))
          dimer = make_shared<Dimer>(itree, ref);
        else
          throw runtime_error("dimerize needs a reference calculation (for now)");
      } else if (form == "r" || form == "refs") {
        vector<shared_ptr<const Reference>> dimer_refs;
        auto units = itree->get_vector<string>("refs", 2);
        for (auto& ikey : units) {
          auto tmp = saved.find(ikey);
          if (tmp == saved.end()) throw runtime_error(string("No reference found with name: ") + ikey);
          else dimer_refs.push_back(static_pointer_cast<const Reference>(tmp->second));
        }

        dimer = make_shared<Dimer>(itree, dimer_refs.at(0), dimer_refs.at(1));
      } else if (form == "linked") {
        dimer = make_shared<Dimer>(itree, ref, /*linked=*/true);
      }

      dimer->scf(itree);

      geom = dimer->sgeom();
      ref = dimer->sref();
    } else if (title == "asd") {
      auto asd = construct_ASD(itree, dimer);
      asd->compute();
    } else if (title == "asd_orbitaloptimize" || title == "asd_orbopt") {
      auto asd = construct_ASD_OrbOpt(itree, dimer);
      asd->compute();
      ref = dimer->sref();
    } else if (title == "multisite") {
      auto multisite = make_shared<MultiSite>(itree, ref);
      ref = multisite->mref();
    } else if (title == "asd_dmrg") {
      auto asd_dmrg = make_shared<RASD>(itree, ref);
      asd_dmrg->project_active();
      asd_dmrg->sweep();
      ref = asd_dmrg->sref();
    } else if (title == "localize") {
      if (ref == nullptr) throw runtime_error("Localize needs a reference");
      if (ref->coeffB()) throw runtime_error("Localize is not implemented for UHF/ROHF");

      string localizemethod = itree->get<string>("algorithm", "pm");
      shared_ptr<OrbitalLocalization> localization;
      if (localizemethod == "region") {
        localization = make_shared<RegionLocalization>(itree, ref);
      }
      else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
        localization = make_shared<PMLocalization>(itree, ref);
      else throw runtime_error("Unrecognized orbital localization method");

      shared_ptr<const Coeff> new_coeff = make_shared<const Coeff>(*localization->localize());
      ref = make_shared<const Reference>(*ref, new_coeff);

    } else if (title == "print") {

      const bool orbitals = itree->get<bool>("orbitals", false);
      const bool vibration = itree->get<bool>("vibration", false);
      const string out_file = itree->get<string>("file", "out.molden");

      if (mpi__->rank() == 0) {
        MoldenOut mfs(out_file);
        mfs << geom;
        if (orbitals || vibration) mfs << ref;
      }
    } else {
      // otherwise, they are considered single point energy calculation
      tie(ignore, ref) = get_energy(title, itree, geom, ref);
    }

    // Save functionality
    string saveref = itree->get<string>("saveref", "");
    if (saveref != "") { saved.insert(make_pair(saveref, static_pointer_cast<const void>(ref))); }

    cout << endl;
    mpi__->barrier();
    timer.tick_print("Method: " + title);
    cout << endl;

  }

  print_footer();
}


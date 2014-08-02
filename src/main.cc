//
// BAGEL - Parallel electron correlation program.
// Filename: main.cc
// Copyright (C) 2009 Toru Shiozaki
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

#include <src/global.h>
#include <src/io/moldenout.h>
#include <src/mp2/mp2grad.h>
#include <src/opt/optimize.h>
#include <src/london/reference_london.h>
#include <src/molecule/localization.h>
#include <src/asd/construct_asd.h>
#include <src/asd_dmrg/rasd.h>
#include <src/util/archive.h>

// debugging
extern void test_solvers(std::shared_ptr<bagel::Geometry>);
extern void test_mp2f12();

// function (TODO will be moved to an appropriate place)

using namespace std;
using namespace bagel;

int main(int argc, char** argv) {

  static_variables();
  print_header();

  try {

    const bool input_provided = argc == 2;
    if (!input_provided) {
      throw runtime_error("no input file provided");
    }
    const string input = argv[1];

    auto idata = make_shared<const PTree>(input);

    shared_ptr<Method> method;
    shared_ptr<Geometry> geom;
    shared_ptr<const Reference> ref;
    shared_ptr<Dimer> dimer;

    map<string, shared_ptr<const void>> saved;
    bool dodf = true;

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
      } else {
        if (!geom) throw runtime_error("molecule block is missing");
        if (!itree->get<bool>("df",true)) dodf = false;
        if (dodf && !geom->df()) throw runtime_error("It seems that DF basis was not specified in molecule block");
      }

      if ((title == "smith" || title == "fci") && ref == nullptr)
        throw runtime_error(title + " needs a reference");

      // most methods are constructed here
      method = construct_method(title, itree, geom, ref);

#ifndef DISABLE_SERIALIZATION
      if (title == "continue") {
        IArchive archive(itree->get<string>("archive"));
        Method* ptr;
        archive >> ptr;
        method = shared_ptr<Method>(ptr);
      }
#endif

      if (method) {

        method->compute();
        ref = method->conv_to_ref();

      } else if (title == "optimize") {

        auto opt = make_shared<Optimize>(itree, geom);
        opt->compute();


      } else if (title == "dimerize") { // dimerize forms the dimer object, does a scf calculation, and then localizes
        const string form = itree->get<string>("form", "displace");
        if (form == "d" || form == "disp" || form == "displace") {
          if (static_cast<bool>(ref))
            dimer = make_shared<Dimer>(itree, ref);
          else
            throw runtime_error("dimerize needs a reference calculation (for now)");
        }
        else if (form == "r" || form == "refs") {
          vector<shared_ptr<const Reference>> dimer_refs;
          auto units = itree->get_vector<string>("refs", 2);
          for (auto& ikey : units) {
            auto tmp = saved.find(ikey);
            if (tmp == saved.end()) throw runtime_error(string("No reference found with name: ") + ikey);
            else dimer_refs.push_back(static_pointer_cast<const Reference>(tmp->second));
          }

          dimer = make_shared<Dimer>(itree, dimer_refs.at(0), dimer_refs.at(1));
        }

        dimer->scf(itree);

        *geom = *dimer->sgeom();
        ref = dimer->sref();
      } else if (title == "asd") {
          auto asd = construct_ASD(itree, dimer);
          asd->compute();
      } else if (title == "asd_dmrg") {
          auto asd = make_shared<RASD>(itree, dimer);
          asd->compute();
      } else if (title == "localize") {
        if (ref == nullptr) throw runtime_error("Localize needs a reference");

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
        const string out_file = itree->get<string>("file", "out.molden");

        if (mpi__->rank() == 0) {
          MoldenOut mfs(out_file);
          mfs << geom;
          if (orbitals) mfs << ref;
        }

      } else {
        if (title != "molecule")
          throw runtime_error("unknown method");
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

  } catch (const exception &e) {
    cout << "  ERROR: EXCEPTION RAISED:" << e.what() << endl;
    throw;
  } catch (...) {
    throw;
  }

  return 0;
}


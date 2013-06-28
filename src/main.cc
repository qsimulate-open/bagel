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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <vector>
#include <tuple>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <memory>

#include <src/scf/overlap.h>
#include <src/scf/coeff.h>
#include <src/wfn/geometry.h>
#include <src/dimer/dimer.h>
#include <src/dimer/dimer_cispace.h>
#include <src/scf/rohf.h>
#include <src/ks/ks.h>
#include <src/io/moldenout.h>
#include <src/wfn/reference.h>
#include <src/rel/relreference.h>
#include <src/wfn/ciwfn.h>
#include <src/fci/distfci.h>
#include <src/fci/harrison.h>
#include <src/fci/knowles.h>
#include <src/casscf/superci.h>
#include <src/casscf/werner.h>
#include <src/casscf/casbfgs.h>
#include <src/mp2/mp2grad.h>
#include <src/global.h>
#include <src/opt/optimize.h>
#include <src/util/constants.h>
#include <src/util/localization.h>
#include <src/util/timer.h>
#include <src/util/lexical_cast.h>
#include <src/rel/dirac.h>
#include <src/rel/relfci.h>
#include <src/transp/transp.h>
#include <src/smith/smith.h>
#include <src/meh/meh.h>

#include <config.h>

// input parser
#include <src/util/input.h>


// debugging
extern void test_solvers(std::shared_ptr<bagel::Geometry>);
extern void test_mp2f12();

using namespace std;
using namespace bagel;

int main(int argc, char** argv) {

  static_variables(argc, argv);
  print_header();

  try {

    const bool input_provided = argc == 2;
    if (!input_provided) {
      throw runtime_error("no input file provided");
    }
    const string input = argv[1];

    auto idata = make_shared<const PTree>(input);

    bool scf_done = false;
    bool casscf_done = false;
    shared_ptr<Geometry> geom;
    shared_ptr<SCF_base> scf;
    shared_ptr<const Reference> ref;
    shared_ptr<const RelReference> relref;
    shared_ptr<Dimer> dimer;

    map<string, shared_ptr<const void>> saved;

    // timer for each method
    Timer timer(-1);

    shared_ptr<const PTree> keys = idata->get_child("bagel");

    for (auto& itree : *keys) {

      string method = itree->get<string>("title", "");
      transform(method.begin(), method.end(), method.begin(), ::tolower);
      if (method.empty()) throw runtime_error("title is missing in one of the input blocks");

      if (method == "molecule") {
        if (ref != nullptr) geom->discard_df();
        geom = make_shared<Geometry>(itree);
        if (itree->get<bool>("restart", false)) {
          ref = shared_ptr<const Reference>();
          relref = shared_ptr<const RelReference>();
        }
        if (ref != nullptr) ref = ref->project_coeff(geom);
        if (relref != nullptr) relref = relref->project_coeff(geom);
      } else {
        if (geom == nullptr) throw runtime_error("molecule block is missing");
      }

      // TODO we are checking this for non-DF methods...
      if (geom->df() == nullptr)
        throw runtime_error("It seems that DF basis was not specified in Geometry");

      if (method == "hf") {

        scf = make_shared<SCF>(itree, geom, ref);
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "dhf") {

        auto dirac = relref ? make_shared<Dirac>(itree, geom, relref) : make_shared<Dirac>(itree, geom, ref);
        dirac->compute();
        relref = dirac->conv_to_ref();

      } else if (method == "relfci") {
        //currently under construction
        auto dirac = relref ? make_shared<Dirac>(itree, geom, relref) : make_shared<Dirac>(itree, geom, ref);
        dirac->compute();
        relref = dirac->conv_to_ref();

        auto relfci = make_shared<RelFCI>(itree, geom, relref);
        relfci->compute();

      } else if (method == "ks") {

        scf = make_shared<KS>(itree, geom, ref);
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "uhf") {

        scf = make_shared<UHF>(itree, geom, ref);
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "rohf") {

        scf = make_shared<ROHF>(itree, geom, ref);
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "optimize") {

        auto opt = make_shared<Optimize>(itree, geom);
        opt->compute();

      } else if (method == "casscf") {

        shared_ptr<CASSCF> casscf;
        string algorithm = itree->get<string>("algorithm", "");
        if (algorithm == "superci" || algorithm == "") {
          casscf = make_shared<SuperCI>(itree, geom, ref);
        } else if (algorithm == "werner" || algorithm == "knowles") {
          casscf = make_shared<WernerKnowles>(itree, geom);
        } else if (algorithm == "bfgs") {
          casscf = make_shared<CASBFGS>(itree, geom, ref);
        } else {
          stringstream ss; ss << "unknown CASSCF algorithm specified: " << algorithm;
          throw runtime_error(ss.str());
        }
        casscf->compute();
        ref = casscf->conv_to_ref();

      } else if (method == "mp2") {

        auto mp2 = make_shared<MP2>(itree, geom);
        mp2->compute();

      } else if (method == "transp") {

        auto tran = make_shared<Transp>(itree, geom, ref);
        tran->compute();

      } else if (method == "smith") {

        if (ref == nullptr) throw runtime_error("SMITH needs a reference");
        auto smith = make_shared<Smith>(itree, ref);
        smith->compute();

      } else if (method == "fci") {
        if (ref == nullptr) throw runtime_error("FCI needs a reference");
        shared_ptr<FCI> fci;

        string algorithm = itree->get<string>("algorithm", "");
        if (algorithm == "" || algorithm == "auto") {
          // TODO At the moment this doesn't take freezing of orbitals into account
          const int nele = ref->geom()->nele();
          const int norb = ref->geom()->nbasis();
          if ( nele <= norb ) fci = make_shared<HarrisonZarrabian>(itree, ref);
          else fci = make_shared<KnowlesHandy>(itree, ref);
        } else if (algorithm == "kh" || algorithm == "knowles" || algorithm == "handy") {
          fci = make_shared<KnowlesHandy>(itree, ref);
        } else if (algorithm == "hz" || algorithm == "harrison" || algorithm == "zarrabian") {
          fci = make_shared<HarrisonZarrabian>(itree, ref);
#ifdef HAVE_MPI_H
        } else if (algorithm == "parallel" || algorithm == "dist") {
          fci = make_shared<DistFCI>(itree, ref);
#endif
        } else {
          throw runtime_error("unknown FCI algorithm specified.");
        }

        fci->compute();

      } else if (method == "dimerize") { // dimerize forms the dimer object, does a scf calculation, and then localizes
        shared_ptr<const PTree> dimdata = itree;

        const string form = dimdata->get<string>("form", "displace");
        if (form == "d" || form == "disp" || form == "displace") {
          double scale = (dimdata->get<bool>("angstrom", false) ? ang2bohr__ : 1.0 ) ;

          double dx = dimdata->get<double>("dx", 0.0) * scale;
          double dy = dimdata->get<double>("dy", 0.0) * scale;
          double dz = dimdata->get<double>("dz", 0.0) * scale;
          array<double,3> disp = {{dx,dy,dz}};

          if (static_cast<bool>(ref)) {
            dimer = make_shared<Dimer>(ref, disp);
          }
          else {
            throw runtime_error("dimerize needs a reference calculation (for now)");
          }
        }
        else if (form == "r" || form == "refs") {
          vector<shared_ptr<const Reference>> dimer_refs;
          auto units = itree->get_child("refs");
          if (units->size() != 2) throw runtime_error("Must provide exactly two references to dimerize with references");
          for (auto& i : *units) {
            string istring = lexical_cast<string>(i->data());
            auto tmp = saved.find(istring);
            if (tmp == saved.end()) throw runtime_error(string("No reference found with name: ") + istring);
            else dimer_refs.push_back(static_pointer_cast<const Reference>(tmp->second));
          }

          dimer = make_shared<Dimer>(dimer_refs.at(0), dimer_refs.at(1));
        }

        dimer->scf(itree);

        *geom = *dimer->sgeom();
        ref = dimer->sref();
      } else if (method == "to-dimer") {
#if 0
        auto inputdata = iter->second;
        auto stringiter = inputdata.find("dimer_active");
        if (stringiter == inputdata.end()) throw runtime_error("Need to set dimer_active in to-dimer method");
        ref = ref->set_active(stringiter->second);

        vector<int> sizes;
        auto bound = iter->second.equal_range("unit");
        for (auto isizes = bound.first; isizes != bound.second; ++isizes) sizes.push_back(lexical_cast<int>(isizes->second));

        dimer = make_shared<Dimer>(ref, make_pair(sizes.at(0), sizes.at(1)));
#else
throw logic_error("broken!");
#endif
      } else if (method == "meh") {
          shared_ptr<DimerCISpace> cispace = dimer->compute_cispace(itree);

          auto meh = make_shared<MultiExcitonHamiltonian>(itree, dimer, cispace);
          meh->compute();
      } else if (method == "localize") {
        if (ref == nullptr) throw runtime_error("Localize needs a reference");

        string localizemethod = itree->get<string>("algorithm", "pm");
        shared_ptr<OrbitalLocalization> localization;
        if (localizemethod == "region") {
          vector<int> sizes;
          auto sizedata = itree->get_child("region_sizes");
          if (sizedata) {
            for (auto& isize : *sizedata) sizes.push_back(lexical_cast<double>(isize->data()));
          }
          else throw runtime_error("Must specify region_sizes with region localization");

          localization = make_shared<RegionLocalization>(ref, sizes);
        }
        else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
          localization = make_shared<PMLocalization>(ref);
        else throw runtime_error("Unrecognized orbital localization method");

        const int max_iter = itree->get<int>("max_iter", 50);
        const double thresh = itree->get<double>("thresh", 1.0e-6);

        shared_ptr<const Coeff> new_coeff = make_shared<const Coeff>(*localization->localize(max_iter,thresh));
        ref = make_shared<const Reference>(ref, new_coeff);

      } else if (method == "print") {

        const bool orbitals = itree->get<bool>("orbitals", false);
        const string out_file = itree->get<string>("file", "out.molden");

        MoldenOut mfs(out_file);
        mfs << geom;
        if(orbitals) mfs << ref;
        mfs.close();

      }

      // Save functionality
      string saveref = itree->get<string>("saveref", "");
      if (saveref != "") { saved.insert(make_pair(saveref, static_pointer_cast<const void>(ref))); }

      cout << endl;
      timer.tick_print("Method: " + method);
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


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
#include <src/parallel/resources.h>
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
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <src/parallel/mpi_interface.h>

#include <config.h>

// input parser
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>


// TODO they are ugly
// TODO to be determined by the number of threads passed by the arguments --num_threads=8 ?
namespace bagel {
  Resources* resources__;
  MPI_Interface* mpi__;
}

// debugging
extern void test_solvers(std::shared_ptr<bagel::Geometry>);
extern void test_mp2f12();

using namespace std;
using namespace bagel;

int main(int argc, char** argv) {
  // setup MPI interface. It does nothing for serial runs
  mpi__ = new MPI_Interface(argc, argv);
  {
    string snum_threads = getenv_multiple("BAGEL_NUM_THREADS", "OMP_NUM_THREADS");
    const int num_threads = snum_threads.empty() ? thread::hardware_concurrency() : lexical_cast<int>(snum_threads);
    if (num_threads < 1)
      throw runtime_error("Set BAGEL_NUM_THREADS for the number of threads used");
    if (mpi__->rank() == 0)
      cout << "  * using " << num_threads << " threads per process" << endl;
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
    resources__ = new Resources(num_threads);
  }

  print_header();

  try {

    const bool input_provided = argc == 2;
    if (!input_provided) {
      throw runtime_error("no input file provided");
    }
    const string input = argv[1];

    boost::property_tree::ptree idata;
    boost::property_tree::json_parser::read_json(input, idata);

    bool scf_done = false;
    bool casscf_done = false;
    shared_ptr<Geometry> geom;
    multimap<string, shared_ptr<Geometry>> saved_geos;
    shared_ptr<SCF_base> scf;
    shared_ptr<const Reference> ref;
    multimap<string, shared_ptr<const Reference>> saved_refs;
    shared_ptr<const RelReference> relref;
    shared_ptr<Dimer> dimer;

    // timer for each method
    Timer timer(-1);

    boost::property_tree::ptree keys = idata.get_child("bagel"); 
    for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
      string method = iter->second.get<string>("title", "");
      transform(method.begin(), method.end(), method.begin(), ::tolower);
      if (method.empty()) throw runtime_error("title is missing in one of the input blocks");

      if (method == "molecule") {
        if (ref != nullptr) geom->discard_df(); 
        geom = make_shared<Geometry>(iter->second);
        if (iter->second.get<bool>("restart", false)) {
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

        scf = make_shared<SCF>(iter->second, geom, ref);
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "dhf") {

        auto dirac = relref ? make_shared<Dirac>(iter->second, geom, relref)
                            : make_shared<Dirac>(iter->second, geom, ref);
        dirac->compute();
        relref = dirac->conv_to_ref();

      } else if (method == "relfci") {
        //currently under construction
        auto dirac = relref ? make_shared<Dirac>(iter->second, geom, relref)
                            : make_shared<Dirac>(iter->second, geom, ref);
        dirac->compute();
        relref = dirac->conv_to_ref();

        auto relfci = make_shared<RelFCI>(iter->second, geom, relref);
        relfci->compute();

      } else if (method == "ks") {

        scf = make_shared<KS>(iter->second, geom, ref);
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "uhf") {

        scf = make_shared<UHF>(iter->second, geom, ref);
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "rohf") {

        scf = make_shared<ROHF>(iter->second, geom, ref);
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "optimize") {

        auto opt = make_shared<Optimize>(iter->second, geom);
        opt->compute();

      } else if (method == "casscf") {

        shared_ptr<CASSCF> casscf;
        string algorithm = iter->second.get<string>("algorithm", "");
        if (algorithm == "superci" || algorithm == "") {
          casscf = make_shared<SuperCI>(iter->second, geom, ref);
        } else if (algorithm == "werner" || algorithm == "knowles") {
          casscf = make_shared<WernerKnowles>(iter->second, geom);
        } else if (algorithm == "bfgs") {
          casscf = make_shared<CASBFGS>(iter->second, geom, ref);
        } else {
          stringstream ss; ss << "unknown CASSCF algorithm specified: " << algorithm; 
          throw runtime_error(ss.str());
        }
        casscf->compute();
        ref = casscf->conv_to_ref();

      } else if (method == "mp2") {

        auto mp2 = make_shared<MP2>(iter->second, geom);
        mp2->compute();

      } else if (method == "transp") {

        auto tran = make_shared<Transp>(iter->second, geom, ref);
        tran->compute();

      } else if (method == "smith") {

        if (ref == nullptr) throw runtime_error("SMITH needs a reference");
        auto smith = make_shared<Smith>(iter->second, ref);
        smith->compute();

      } else if (method == "fci") {
        if (ref == nullptr) throw runtime_error("FCI needs a reference");
        shared_ptr<FCI> fci;

        string algorithm = iter->second.get<string>("algorithm", "");
        if (algorithm == "" || algorithm == "auto") {
          // TODO At the moment this doesn't take freezing of orbitals into account
          const int nele = ref->geom()->nele();
          const int norb = ref->geom()->nbasis();
          if ( nele <= norb ) fci = make_shared<HarrisonZarrabian>(iter->second, ref);
          else fci = make_shared<KnowlesHandy>(iter->second, ref);
        } else if (algorithm == "kh" || algorithm == "knowles" || algorithm == "handy") {
          fci = make_shared<KnowlesHandy>(iter->second, ref);
        } else if (algorithm == "hz" || algorithm == "harrison" || algorithm == "zarrabian") {
          fci = make_shared<HarrisonZarrabian>(iter->second, ref);
#ifdef HAVE_MPI_H
        } else if (algorithm == "parallel" || algorithm == "dist") {
          fci = make_shared<DistFCI>(iter->second, ref);
#endif
        } else {
          throw runtime_error("unknown FCI algorithm specified.");
        }

        fci->compute();

      } else if (method == "dimerize") { // dimerize forms the dimer object, does a scf calculation, and then localizes
        const boost::property_tree::ptree dimdata = iter->second;

        const string form = dimdata.get<string>("form", "displace");
        if (form == "d" || form == "disp" || form == "displace") {
          double scale = (dimdata.get<bool>("angstrom", false) ? ang2bohr__ : 1.0 ) ;

          double dx = dimdata.get<double>("dx", 0.0) * scale;
          double dy = dimdata.get<double>("dy", 0.0) * scale;
          double dz = dimdata.get<double>("dz", 0.0) * scale;
          array<double,3> disp = {{dx,dy,dz}};

          if (static_cast<bool>(ref)) {
            dimer = make_shared<Dimer>(ref, disp);
          }
          else {
            throw runtime_error("dimerize needs a reference calculation (for now)");
          }
        }
        else if (form == "r" || form == "refs") {
#if 0
          shared_ptr<const Reference> refA;
          shared_ptr<const Reference> refB;

          auto iterA = dimdata.find("unita");
          if ( iterA != dimdata.end() ) {
            auto irefA = saved_refs.find(iterA->second);
            if ( irefA != saved_refs.end() ) {
              refA = irefA->second;
            }
            else {
              throw runtime_error("Saved reference \"" + iterA->second + "\" not found");
            }
          }
          else {
            throw runtime_error("No input provided for unit A");
          }

          auto iterB = dimdata.find("unitb");
          if ( iterB != dimdata.end() ) {
            auto irefB = saved_refs.find(iterB->second);
            if ( irefB != saved_refs.end() ) {
              refB = irefB->second;
            }
            else {
              throw runtime_error("Saved reference \"" + iterB->second + "\" not found");
            }
          }
          else {
            throw runtime_error("No input provided for unit B");
          }

          dimer = make_shared<Dimer>(refA, refB);
#else
throw logic_error("broken!");
#endif
        }

        dimer->scf(iter->second);

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
          shared_ptr<DimerCISpace> cispace = dimer->compute_cispace(iter->second);
    
          auto meh = make_shared<MultiExcitonHamiltonian>(iter->second, dimer, cispace);
          meh->compute();
          meh->print();
      } else if (method == "localize") {
        if (ref == nullptr) throw runtime_error("Localize needs a reference");

        string localizemethod = iter->second.get<string>("algorithm", "pm");
        shared_ptr<OrbitalLocalization> localization;
        if (localizemethod == "region") {
#if 0
          vector<int> sizes;
          auto bound = iter->second.equal_range("region");
          for (auto isizes = bound.first; isizes != bound.second; ++isizes) sizes.push_back(lexical_cast<int>(isizes->second));

          localization = make_shared<RegionLocalization>(ref, sizes);
#else
throw logic_error("broken!");
#endif
        }
        else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
          localization = make_shared<PMLocalization>(ref);
        else throw runtime_error("Unrecognized orbital localization method");

        const int max_iter = iter->second.get<int>("max_iter", 50);
        const double thresh = iter->second.get<double>("thresh", 1.0e-6);

        shared_ptr<const Coeff> new_coeff = make_shared<const Coeff>(*localization->localize(max_iter,thresh));
        ref = make_shared<const Reference>(ref, new_coeff);
        
      } else if (method == "print") {

        const boost::property_tree::ptree pdata = iter->second;
        const bool orbitals = pdata.get<bool>("orbitals", false);
        const string out_file = pdata.get<string>("file", "out.molden");

        MoldenOut mfs(out_file);
        mfs << geom;
        if(orbitals) mfs << ref;
        mfs.close();

      } else if (method == "save") {
#if 0
        auto sdata = iter->second;

        auto igeom = sdata.find("geom");
          if ( igeom != sdata.end() ) saved_geos.insert(make_pair(igeom->second, geom));

        auto iref = sdata.find("ref");
          if ( iref != sdata.end() ) saved_refs.insert(make_pair(iref->second, ref));
#else
throw logic_error("broken!");
#endif
      }
      #if 0 // <---- Testing environment
      else if (method == "testing") {
        multimap<string, string> geominfo = idata->get_input("molecule");
        multimap<string,string> dimdata = iter->second;

        shared_ptr<FCI> fci = make_shared<HarrisonZarrabian>(iter->second, ref);
        fci->compute();

        auto ciwfn = fci->conv_to_ciwfn();

        double scale = (read_input<bool>(dimdata,"angstrom",false) ? ang2bohr__ : 1.0 ) ;

        double dx = read_input<double>(dimdata,"dx",0.0) * scale;
        double dy = read_input<double>(dimdata,"dy",0.0) * scale;
        double dz = read_input<double>(dimdata,"dz",0.0) * scale;
        array<double,3> disp = {{dx,dy,dz}};

        dimer = make_shared<Dimer>(ciwfn, disp);

      }
      #endif

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
  delete resources__;
  delete mpi__;


  return 0;
}


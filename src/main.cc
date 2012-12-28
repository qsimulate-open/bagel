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
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <memory>

#include <boost/lexical_cast.hpp>

#include <src/scf/overlap.h>
#include <src/scf/coeff.h>
#include <src/scf/geometry.h>
#include <src/dimer/dimer.h>
#include <src/dimer/dimer_scf.h>
#include <src/scf/rohf.h>
#include <src/io/moldenout.h>
#include <src/wfn/reference.h>
#include <src/wfn/ciwfn.h>
#include <src/fci/harrison.h>
#include <src/fci/knowles.h>
#include <src/casscf/superci.h>
#include <src/casscf/werner.h>
#include <src/mp2/mp2grad.h>
#include <src/global.h>
#include <src/parallel/resources.h>
#include <src/opt/opt.h>
#include <src/util/input.h>
#include <src/util/constants.h>
#include <src/util/localization.h>
#include <src/rel/dirac.h>
#include <src/smith/storage.h>
#include <src/smith/MP2.h>
#include <src/smith/CAS_all_active.h>
#include <src/smith/CAS_test.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <src/parallel/mpi_interface.h>

#include <config.h>


// TODO they are ugly
// TODO to be determined by the number of threads passed by the arguments --num_threads=8 ?
namespace bagel {
  Resources* resources__;
  MPI_Interface* mpi__;
}

// debugging
extern void test_solvers(std::shared_ptr<bagel::Geometry>);
extern void test_mp2f12();

using std::cout;
using std::endl;
using namespace bagel;

#include <src/parallel/paramatrix.h>

int main(int argc, char** argv) {

  // setup MPI interface. It does nothing for serial runs
  mpi__ = new MPI_Interface(argc, argv);
  {
    // TODO will be interfaced to input
    int num_threads = 16;
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif
    resources__ = new Resources(num_threads);
  }

  print_header();

  try {

    const bool input_provided = argc == 2;
    if (!input_provided) {
      throw std::runtime_error("no input file provided");
    }
    const std::string input = argv[1];

    std::shared_ptr<InputData> idata(new InputData(input));

    bool scf_done = false;
    bool casscf_done = false;
    std::shared_ptr<Geometry> geom;
    std::shared_ptr<SCF_base> scf;
    std::shared_ptr<const Reference> ref;
    std::shared_ptr<Dimer> dimer;

    std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

    for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
      const std::string method = iter->first;

      if (method == "molecule") {
        geom = std::shared_ptr<Geometry>(); // kill the previous one first
        geom = std::shared_ptr<Geometry>(new Geometry(iter->second));
        if (read_input<bool>(iter->second, "restart", false)) ref = std::shared_ptr<const Reference>();
        if (ref != nullptr) ref = ref->project_coeff(geom);
      } else {
        if (geom == nullptr) throw std::runtime_error("molecule block is missing");
      }

      if (method.substr(0,3) == "df-" && geom->df() == nullptr)
        throw std::runtime_error("It seems that DF basis was not specified in Geometry");

      if (method == "hf") {

        scf = std::shared_ptr<SCF<0> >(new SCF<0>(iter->second, geom, ref));
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "dhf") {

        scf = std::shared_ptr<Dirac>(new Dirac(iter->second, geom, ref));
        scf->compute();
//      ref = scf->conv_to_ref();

      } else if (method == "df-hf") {

        scf = std::shared_ptr<SCF<1> >(new SCF<1>(iter->second, geom, ref));
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "df-uhf" || method == "uhf") {

        scf = std::shared_ptr<UHF>(new UHF(iter->second, geom, ref));
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "df-rohf" || method == "rohf") {

        scf = std::shared_ptr<ROHF>(new ROHF(iter->second, geom, ref));
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "df-uhf-opt" || method == "uhf-opt") {

        std::shared_ptr<Opt<UHF> > opt(new Opt<UHF>(idata, iter->second, geom));
        for (int i = 0; i != 100; ++i)
          if (opt->next()) break;

      } else if (method == "df-rohf-opt" || method == "rohf-opt") {

        std::shared_ptr<Opt<ROHF> > opt(new Opt<ROHF>(idata, iter->second, geom));
        for (int i = 0; i != 100; ++i)
          if (opt->next()) break;

      } else if (method == "df-hf-opt") {

        std::shared_ptr<Opt<SCF<1> > > opt(new Opt<SCF<1> >(idata, iter->second, geom));
        for (int i = 0; i != 100; ++i)
          if (opt->next()) break;

      } else if (method == "casscf") {

        std::shared_ptr<CASSCF> casscf;
        std::string algorithm = read_input<std::string>(iter->second, "algorithm", "");
        if (algorithm == "superci" || algorithm == "") {
          casscf = std::shared_ptr<CASSCF>(new SuperCI(iter->second, geom));
        } else if (algorithm == "werner" || algorithm == "knowles") {
          casscf = std::shared_ptr<CASSCF>(new WernerKnowles(iter->second, geom));
        } else {
          throw std::runtime_error("unknown CASSCF algorithm specified.");
        }
        casscf->compute();
        ref = casscf->conv_to_ref();

      } else if (method == "casscf-opt") {
        std::string algorithm = read_input<std::string>(iter->second, "algorithm", "");
        // in case of SS-CASSCF
        if (read_input<int>(iter->second, "nstate", 1) == 1) {
          if (algorithm == "superci" || algorithm == "") {
            std::shared_ptr<Opt<SuperCI> > opt(new Opt<SuperCI>(idata, iter->second, geom));
            for (int i = 0; i != 100; ++i)
              if (opt->next()) break;
          } else if (algorithm == "werner" || algorithm == "knowles") {
            std::shared_ptr<Opt<WernerKnowles> > opt(new Opt<WernerKnowles>(idata, iter->second, geom));
            for (int i = 0; i != 100; ++i)
              if (opt->next()) break;
          } else {
            throw std::runtime_error("unknown CASSCF algorithm specified.");
          }
        // in case of SA-CASSCF
        } else {
          if (algorithm == "superci" || algorithm == "") {
            std::shared_ptr<Opt<SuperCIGrad> > opt(new Opt<SuperCIGrad>(idata, iter->second, geom));
            for (int i = 0; i != 100; ++i)
              if (opt->next()) break;
          } else {
            throw std::runtime_error("unknown CASSCF algorithm specified.");
          }


        }
      } else if (method == "mp2") {

        std::shared_ptr<MP2> mp2(new MP2(iter->second, geom));
        mp2->compute();

      } else if (method == "mp2-opt") {

        std::shared_ptr<Opt<MP2Grad> > opt(new Opt<MP2Grad>(idata, iter->second, geom));
        for (int i = 0; i != 100; ++i)
          if (opt->next()) break;

      } else if (method == "smith") {

        std::string method = read_input<std::string>(iter->second, "method", "mp2");
        if (ref == nullptr) throw std::runtime_error("SMITH needs a reference");
        if (method == "mp2") {
          std::shared_ptr<SMITH::MP2::MP2<SMITH::Storage_Incore> > mp2(new SMITH::MP2::MP2<SMITH::Storage_Incore>(ref));
          mp2->solve();
        } else if (method == "caspt2") {
          std::shared_ptr<SMITH::CAS_all_active::CAS_all_active<SMITH::Storage_Incore> > cas(new SMITH::CAS_all_active::CAS_all_active<SMITH::Storage_Incore>(ref));
          cas->solve();
        } else if (method == "caspt2-test") {
          std::shared_ptr<SMITH::CAS_test::CAS_test<SMITH::Storage_Incore> > cas(new SMITH::CAS_test::CAS_test<SMITH::Storage_Incore>(ref));
          cas->solve();
        } else {
          std::stringstream ss; ss << method << " method is not implemented in SMITH";
          throw std::logic_error(ss.str());
        }

      } else if (method == "fci") {
        if (ref == nullptr) throw std::runtime_error("FCI needs a reference");
        std::shared_ptr<FCI> fci;

        std::string algorithm = read_input<std::string>(iter->second, "algorithm", "");
        if (algorithm == "" || algorithm == "auto") {
          // TODO At the moment this doesn't take freezing of orbitals into account
          const int nele = ref->geom()->nele();
          const int norb = ref->geom()->nbasis();
          if ( (nele) <= norb ) fci = std::shared_ptr<FCI>(new HarrisonZarrabian(iter->second, ref));
          else fci = std::shared_ptr<FCI>(new KnowlesHandy(iter->second, ref));
        } else if (algorithm == "kh" || algorithm == "knowles" || algorithm == "handy") {
          fci = std::shared_ptr<FCI>(new KnowlesHandy(iter->second, ref));
        } else if (algorithm == "hz" || algorithm == "harrison" || algorithm == "zarrabian") {
          fci = std::shared_ptr<FCI>(new HarrisonZarrabian(iter->second, ref));
        } else {
          throw std::runtime_error("unknown FCI algorithm specified.");
        }

        fci->compute();

      } else if (method == "dimerize") {

        std::multimap<std::string,std::string> dimdata = iter->second;

        double scale = (read_input<bool>(dimdata,"angstrom",false) ? ang2bohr__ : 1.0 ) ;

        double dx = read_input<double>(dimdata,"dx",0.0) * scale;
        double dy = read_input<double>(dimdata,"dy",0.0) * scale;
        double dz = read_input<double>(dimdata,"dz",0.0) * scale;
        std::array<double,3> disp = {{dx,dy,dz}};

        if (static_cast<bool>(ref)) {
          dimer = std::shared_ptr<Dimer>(new Dimer(ref,disp));
        }
        else {
          dimer = std::shared_ptr<Dimer>(new Dimer(geom,disp));
        }

        geom = dimer->sgeom();
        ref = dimer->sref();

      #if 0
      } else if (method == "dimer-scf") {

        scf = std::shared_ptr<DimerSCF>(new DimerSCF(iter->second, dimer));
        scf->compute();
        ref = scf->conv_to_ref();
        dimer->set_sref(ref);

        geom = dimer->sgeom();
        ref = dimer->sref();
        
      #endif
      } else if (method == "localize") {
        if (ref == nullptr) throw std::runtime_error("Localize needs a reference");

        std::string localizemethod = read_input<std::string>(iter->second,"algorithm", "pm");
        std::shared_ptr<OrbitalLocalization> localization;
        if (localizemethod == "region") {
          std::vector<int> sizes;
          auto bound = iter->second.equal_range("region");
          for (auto isizes = bound.first; isizes != bound.second; ++isizes) sizes.push_back(boost::lexical_cast<int>(isizes->second));

          localization = std::shared_ptr<OrbitalLocalization>(new RegionLocalization(ref, sizes));
        }
        else if (localizemethod == "pm" || localizemethod == "pipek" || localizemethod == "mezey" || localizemethod == "pipek-mezey")
          localization = std::shared_ptr<OrbitalLocalization>(new PMLocalization(ref));
        else throw std::runtime_error("Unrecognized orbital localization method");

        const int max_iter = read_input<int>(iter->second,"max_iter", 50);
        const double thresh = read_input<double>(iter->second,"thresh", 1.0e-6);

        std::shared_ptr<const Coeff> new_coeff = localization->localize(max_iter,thresh);
        ref = std::shared_ptr<const Reference>(new const Reference( ref, new_coeff ));
        
      } else if (method == "print") {

        std::multimap<std::string, std::string> pdata = iter->second;
        bool orbitals = read_input<bool>(pdata, "orbitals", false);
        std::string out_file = read_input<std::string>(pdata, "file", "out.molden");

        MoldenOut mfs(out_file);
        mfs << geom;
        if(orbitals) mfs << ref;
        mfs.close();

      }
      #if 0 // <---- Testing environment
      else if (method == "testing") {
        std::multimap<std::string, std::string> testdata = idata->get_input("testing");
        std::multimap<std::string, std::string> geominfo = idata->get_input("molecule");

        #if 0
        std::shared_ptr<Matrix> density = ref->coeff()->form_density_rhf(ref->nocc());
        // for testing purposes, only using region_size right now which assigns the first "region_size" atoms to the first
        // region and the rest to the second region.
        int region_size = read_input<int>(testdata, "region_size", 0);
        if(region_size <= 0) throw std::runtime_error("region_size should be greater than 0");
        
        std::pair<int, int> pair1(0, region_size-1); std::pair<int,int> pair2(region_size,geom->natom()-1);
        std::vector<std::pair<int, int> > bounds = {pair1, pair2};
        std::shared_ptr<OrbitalLocalization> regionalize(new RegionLocalization(geom, density, bounds));

        std::shared_ptr<Matrix> regional_mos = regionalize->localize();
        //regional_mos->print("Regionalized MOs", geom->nbasis());

        std::shared_ptr<const Coeff> regional_coeff(new const Coeff(*regional_mos));

        #endif

        #if 1
        std::shared_ptr<Matrix> overlap(new Overlap(geom));
        std::shared_ptr<OrbitalLocalization> pml(new PMLocalization(ref->geom(), ref->coeff(), ref->nocc()));
        std::shared_ptr<Matrix> new_matrix = pml->localize();
        std::shared_ptr<const Coeff> new_coeff(new const Coeff(*new_matrix));
        #endif

        ref = std::shared_ptr<Reference>(new Reference(geom, new_coeff, ref->nclosed(), ref->nact(), ref->nvirt()));
        //dimer->set_sref(ref);
        //dimer->fci(testdata);
        //dimer->hamiltonian();
      }
      #endif
    }
    print_footer();

    //test_solvers(geom);

  } catch (const std::exception &e) {
    cout << "  ERROR: EXCEPTION RAISED:" << e.what() << endl;
    throw;
  } catch (...) {
    throw;
  }
  delete resources__;
  delete mpi__;


  return 0;
}


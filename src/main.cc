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
#include <src/util/resources.h>
#include <src/opt/opt.h>
#include <src/util/input.h>
#include <src/util/constants.h>
#include <src/rel/dirac.h>
#include <src/smith/storage.h>
#include <src/smith/MP2.h>
#include <src/smith/CAS_all_active.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <src/parallel/mpi_interface.h>


using namespace bagel;
// TODO to be determined by the number of threads passed by the arguments --num_threads=8 ?
Resources* bagel::resources__;

// debugging
extern void smith_test(std::shared_ptr<Reference>);
extern void test_solvers(std::shared_ptr<Geometry>);
extern void test_mp2f12();
extern void test_grad(std::shared_ptr<Reference>);

using std::cout;
using std::endl;

int main(int argc, char** argv) {

  // testing MPI.
#if 0
  MPI_Interface mpi_;
#endif

  try {
    print_header();

    {
      // TODO will be interfaced to input
      int num_threads = 8;
#ifdef _OPENMP
      omp_set_num_threads(num_threads);
#endif
      resources__ = new Resources(num_threads);
      cout << std::endl << "  " <<  num_threads << " thread" << (num_threads == 1 ? "" : "s") << " will be used" << std::endl << std::endl;
    }

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

      }
        else if (method == "dimerize") {

        std::multimap<std::string,std::string> dimdata = iter->second;

        double scale = (read_input<bool>(dimdata,"angstrom",false) ? ang2bohr__ : 1.0 ) ;

        double dx = read_input<double>(dimdata,"dx",0.0) * scale;
        double dy = read_input<double>(dimdata,"dy",0.0) * scale;
        double dz = read_input<double>(dimdata,"dz",0.0) * scale;
        std::array<double,3> disp = {{dx,dy,dz}};

        std::shared_ptr<Dimer> dim;
        if (static_cast<bool>(ref)) {
          dim = std::shared_ptr<Dimer>(new Dimer(ref,disp));
        }
        else {
          dim = std::shared_ptr<Dimer>(new Dimer(geom,disp));
        }

        geom = dim->sgeom();
        ref = dim->sref();

      } else if (method == "print") {

        std::multimap<std::string, std::string> pdata = iter->second;
        bool orbitals = read_input<bool>(pdata, "orbitals", false);
        std::string out_file = read_input<std::string>(pdata, "file", "out.molden");

        MoldenOut mfs(out_file);
        mfs << geom;
        if(orbitals) mfs << ref;
        mfs.close();

      }
      #if 1 // <---- Testing environment
      else if (method == "testing") {
        std::multimap<std::string, std::string> testdata = idata->get_input("testing");
        std::multimap<std::string, std::string> geominfo = idata->get_input("molecule");

        //std::shared_ptr<FCI> fci(new HarrisonZarrabian(iter->second, ref));
        //fci->compute();
        //std::shared_ptr<const CIWfn> ci = fci->conv_to_ciwfn();

        double dx = read_input<double>(testdata, "dx", 0.0);
        double dy = read_input<double>(testdata, "dy", 0.0);
        double dz = read_input<double>(testdata, "dz", 0.0);
        std::array<double,3> disp = {{dx,dy,dz}};
        std::shared_ptr<Dimer> dim(new Dimer(ref, disp));
        //dim->hamiltonian();

        std::shared_ptr<DimerSCF<1> > dimer_scf(new DimerSCF<1>(testdata, dim));
        dimer_scf->compute();

        #if 0
        std::shared_ptr<Coeff> ovlp = dim->overlap();
        //ovlp->print("overlap", 12);
        double max = 0.0;
        const int n = ovlp->ndim();
        const int m = ovlp->mdim();
        
        for(int i = 0; i < m; ++i) {
          for(int j = 0; j < n; ++j) {
            if (i==j) continue;
            double value = abs(ovlp->element(i,j));
            if (value > max) max = value;
          }
        }
        std::cout << std::setprecision(12) << "max: " << max << std::endl;
        #endif
      
      }
      #endif
    }
    print_footer();

    /////////////////////////////////////
    //smith_test(ref);
    /////////////////////////////////////
    //test_solvers(geom);

    delete resources__;

  } catch (const std::exception &e) {
    cout << "  ERROR: EXCEPTION RAISED:" << e.what() << endl;
    throw;
  } catch (...) {
    throw;
  }

  return 0;
}


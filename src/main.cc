//
// Newint - Parallel electron correlation program.
// Filename: main.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <src/scf/geometry.h>
#include <src/scf/rohf.h>
#include <src/molden/molden.h>
#include <src/wfn/reference.h>
#include <src/fci/fci.h>
#include <src/casscf/superci.h>
#include <src/casscf/werner.h>
#include <src/mp2/mp2grad.h>
#include <src/global.h>
#include <src/stackmem.h>
#include <src/opt/opt.h>
#include <src/util/input.h>

StackMem* stack;

// debugging
extern void smith_test(std::shared_ptr<Reference>);
extern void test_solvers(std::shared_ptr<Geometry>);
extern void test_mp2f12();
extern void test_grad(std::shared_ptr<Reference>);

using std::cout;
using std::endl;

int main(int argc, char** argv) {
  // openmp is broken now due to the use of stack.
  // What we need is a proper thread model.
  #ifdef _OPENMP
  assert(false); // trap
  #endif

#if 0
  stack = new StackMem(static_cast<size_t>(100000000));
  test_mp2f12();
  abort();
#endif

  try {
    print_header();

    const bool input_provided = argc == 2;
    if (!input_provided) {
      throw std::runtime_error("no input file provided");
    }
    const std::string input = argv[1];

    std::shared_ptr<InputData> idata(new InputData(input));

    { // stack
      std::multimap<std::string, std::string> geominfo = idata->get_input("molecule");
      double size = 1.0e6;
      auto iter = geominfo.find("stack");
      if (iter != geominfo.end()) {
        std::string p = iter->second;
        if (p.find("m") != std::string::npos)
          size = 1.0e6*boost::lexical_cast<int>(p.erase(p.size()-1));
        else if (p.find("g") != std::string::npos)
          size = 1.0e9*boost::lexical_cast<int>(p.erase(p.size()-1));
      }
      stack = new StackMem(static_cast<size_t>(size));
      cout << "  Stack memory of " << std::setprecision(2) << std::fixed << size*8.0e-6 << " MB allocated" << endl << endl; 
    }
    std::shared_ptr<Geometry> geom(new Geometry(idata));
    std::list<std::pair<std::string, std::multimap<std::string, std::string> > > keys = idata->data();

    bool scf_done = false;
    bool casscf_done = false;
    std::shared_ptr<SCF_base> scf;
    std::shared_ptr<const Reference> ref;

    for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
      const std::string method = iter->first;

      if (method.substr(0,3) == "df-" && !static_cast<bool>(geom->df()))
        throw std::runtime_error("It seems that DF basis was not specified in Geometry");

      if (method == "hf") {

        scf = std::shared_ptr<SCF<0> >(new SCF<0>(iter->second, geom));
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "df-hf") {

        scf = std::shared_ptr<SCF<1> >(new SCF<1>(iter->second, geom));
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "df-uhf" || method == "uhf") {

        scf = std::shared_ptr<UHF>(new UHF(iter->second, geom));
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "df-rohf" || method == "rohf") {

        scf = std::shared_ptr<ROHF>(new ROHF(iter->second, geom));
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

      } else if (method == "fci") {
        if (!static_cast<bool>(ref)) throw std::runtime_error("FCI needs a reference");

        std::shared_ptr<FCI> fci(new FCI(iter->second, ref));
        fci->compute();

      } else if (method == "mp2") {

        std::shared_ptr<MP2> mp2(new MP2(iter->second, geom));
        mp2->compute();

      } else if (method == "mp2-opt") {

        std::shared_ptr<Opt<MP2Grad> > opt(new Opt<MP2Grad>(idata, iter->second, geom));
        for (int i = 0; i != 100; ++i)
          if (opt->next()) break;

      } else if (method == "print") {
        std::multimap<std::string, std::string> pdata = idata->get_input("print");
        bool orbitals = read_input<bool>(pdata, "orbitals", false);
        std::string out_file = read_input<std::string>(pdata, "file", "out.molden");

        Molden molden;
        molden.write_geo(geom, out_file);
        if (orbitals) molden.write_mos(ref, out_file);
      }
      #if 0 // <---- I've probably messed something up so I'm just taking this out for now
      else if (method == "molden") {

        std::multimap<std::string, std::string> geominfo = idata->get_input("molecule");
        std::string molden_file(idata->get_input("molden").find("molden")->second);

        Molden mf;
        std::shared_ptr<Geometry> molden_geometry(new Geometry( mf.read_geo(molden_file), iter->second) );

        cout << molden_geometry->nuclear_repulsion() << endl;

        scf = std::shared_ptr<SCF<0> >(new SCF<0>(iter->second, molden_geometry));
        scf->compute();
       
      }
      #endif
    }
    print_footer();

    /////////////////////////////////////
    //smith_test(ref);
    /////////////////////////////////////
    //test_solvers(geom);

    delete stack;

  } catch (const std::exception &e) {
    cout << "  ERROR: EXCEPTION RAISED:" << e.what() << endl;
    throw;
  } catch (...) {
    throw;
  }

  return 0;
}


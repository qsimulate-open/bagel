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
#include <src/scf/scf.h>
#include <src/wfn/reference.h>
#include <src/fci/fci.h>
#include <src/casscf/superci.h>
#include <src/casscf/werner.h>
#include <src/mp2/mp2.h>
#include <src/global.h>
#include <src/stackmem.h>
#include <src/util/input.h>

using namespace std;

StackMem* stack;

// debugging
extern void smith_test(shared_ptr<Reference>);
extern void test_solvers(shared_ptr<Geometry>);
extern void test_mp2f12();
extern void test_grad(shared_ptr<Reference>);

int main(int argc, char** argv) {
  // openmp is broken now due to the use of stack.
  // What we need is a proper thread model.
  #ifdef _OPENMP
  assert(false); // trap
  #endif

  //test_mp2f12();
  //abort();

  try {
    print_header();

    const bool input_provided = argc == 2;
    if (!input_provided) {
      throw runtime_error("no input file provided");
    }
    const string input = argv[1];

    shared_ptr<InputData> idata(new InputData(input));

    const bool fci_card = idata->exist("fci"); 
    const bool casscf_card = idata->exist("casscf");

    { // stack
      multimap<string, string> geominfo = idata->get_input("molecule");
      double size = 1.0e6;
      auto iter = geominfo.find("stack");
      if (iter != geominfo.end()) {
        string p = iter->second;
        if (p.find("m") != string::npos)
          size = 1.0e6*boost::lexical_cast<int>(p.erase(p.size()-1));
        else if (p.find("g") != string::npos)
          size = 1.0e9*boost::lexical_cast<int>(p.erase(p.size()-1));
      }
      stack = new StackMem(static_cast<size_t>(size));
      cout << "  Stack memory of " << setprecision(2) << fixed << size*8.0e-6 << " MB allocated" << endl << endl; 
    }
    shared_ptr<Geometry> geom(new Geometry(idata));
    list<pair<string, multimap<string, string> > > keys = idata->data();

    bool scf_done = false;
    bool casscf_done = false;
    shared_ptr<SCF_base> scf;
    shared_ptr<CASSCF> casscf;
    shared_ptr<Reference> ref;

    for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
      const string method = iter->first;

      if (method == "hf") {

        shared_ptr<SCF<0> > scf_(new SCF<0>(iter->second, geom)); scf = scf_;
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "df-hf") {

        if (!geom->df()) throw runtime_error("It seems that DF basis was not specified in Geometry");
        shared_ptr<SCF<1> > scf_(new SCF<1>(iter->second, geom)); scf = scf_;
        scf->compute();
        ref = scf->conv_to_ref();

      } else if (method == "casscf") {
        if (!ref) throw runtime_error("CASSCF needs a reference");

        string algorithm = read_input<string>(iter->second, "algorithm", ""); 
        if (algorithm == "superci" || algorithm == "") {
          shared_ptr<CASSCF> casscf_(new SuperCI(iter->second, geom, ref)); casscf = casscf_;
          casscf->compute();
          ref = casscf->conv_to_ref();
        } else if (algorithm == "werner" || algorithm == "knowles") {
          shared_ptr<CASSCF> werner(new WernerKnowles(iter->second, geom, ref));
          werner->compute();
          ref = werner->conv_to_ref();
        } else {
          throw runtime_error("unknown CASSCF algorithm specified.");
        }

      } else if (method == "fci") {
        if (!ref) throw runtime_error("FCI needs a reference");

        shared_ptr<FCI> fci(new FCI(iter->second, geom, ref));
        fci->compute();

      } else if (method == "mp2") {
        if (!ref) throw runtime_error("MP2 needs a reference");

        shared_ptr<MP2> mp2(new MP2(iter->second, geom, ref));
        mp2->compute();
      }
    }
    print_footer();

    /////////////////////////////////////
    //smith_test(ref);
    /////////////////////////////////////
    //test_solvers(geom);
    test_grad(ref);
    /////////////////////////////////////

    delete stack;

  } catch (const std::exception &e) {
    cout << "  ERROR: EXCEPTION RAISED:" << e.what() << endl;
    throw;
  } catch (...) {
    throw;
  }

  return 0;
}


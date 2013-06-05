//
// BAGEL - Parallel electron correlation program.
// Filename: transp.cc
// Copyright (C) 2013 Toru Shiozaki
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


#include <src/transp/transp.h>
#include <src/scf/scf.h>
#include <src/fci/knowles.h>
#include <src/fci/harrison.h>

using namespace std;
using namespace bagel;

Transp::Transp(const boost::property_tree::ptree& idata_, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> ref)
  : ref_(ref) {

  // for HF, Reference (src/wfn/referenceh.h) has all the information you need.
  ref_->coeff()->print("MO coefficient", 12);

  cout << endl;
  cout << " eigenvalues" << endl;
  int n = 0;
  for (auto& i : ref_->eig()) {
    cout << setw(10) << n++ << setw(20) << setprecision(10) << i << endl;
  }

  // number of state (now default to 15)
  nstate_ = idata_.get<int>("nstate", 15);
}


void Transp::compute() {

  // diagonalization within active space for the ground state
  // 3 orbitals from HOMO, 3 orbitals from LUMO. Unfortunately some pi orbitals are outside this window. One needs to pick by hand (or write a code to pick them).
  {
    KnowlesHandy g(boost::property_tree::ptree(), ref_, 18, 6, 1);
    g.compute();
  }
  // do the same but with 25 states. If you go further it seems numerical noise kills the calculation.
  // Also if you add an input keyword "nstate = 10" within transp, you can change the numer of states
  {
    // all of them are singlet
    boost::property_tree::ptree options;
    options.put("thresh", "1.0e-15");
    KnowlesHandy g(options, ref_, 18, 6, nstate_);
    g.compute();
  }
  {
    // all of them are triplet
    boost::property_tree::ptree options;
    options.put("nspin", "2");
    options.put("thresh", "1.0e-15");
    KnowlesHandy g(options, ref_, 18, 6, nstate_);
    g.compute();
  }
  // you can change the number of electrons
  {
    // ionized states. doublets
    boost::property_tree::ptree options;
    options.put("charge", "1");
    options.put("nspin", "1");
    options.put("thresh", "1.0e-15");
    KnowlesHandy g(options, ref_, 18, 6, nstate_);
    g.compute();
  }
  // coupling term requires a little bit of coding - Shane might already compute it, but it is anyway simple.
}

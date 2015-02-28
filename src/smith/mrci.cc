//
// BAGEL - Parallel electron correlation program.
// Filename: mrci.cc
// Copyright (C) 2015 Toru Shiozaki
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


#include <src/multi/casscf/cashybrid.h>
#include <src/smith/smith.h>
#include <src/smith/mrci.h>


using namespace std;
using namespace bagel;

MRCI::MRCI(shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : Method(inp, geom, ref) {

  Timer timer;

  // compute CASSCF first
  auto cas = make_shared<CASHybrid>(inp, geom, ref);
  cas->compute();

  // update reference
  ref_ = cas->conv_to_ref();
  thresh_ = cas->thresh();
  ref_energy_ = cas->energy();

  timer.tick_print("Reference calculation");

  cout << endl << "  === MRCI calculation ===" << endl << endl;
}


// compute smith and set rdms and ci deriv to a member
void MRCI::compute() {

  // construct SMITH here
  shared_ptr<const PTree> smithinput = idata_->get_child("smith");
  auto smith = make_shared<Smith>(smithinput, ref_->geom(), ref_);
  smith->compute();

}

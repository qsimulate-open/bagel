//
// BAGEL - Parallel electron correlation program.
// Filename: zuperci.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson Bates <jefferson.bates@northwestern.edu>
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


#include <src/zcasscf/zqvec.h>
#include <src/math/step_restrict_bfgs.h>
#include <src/rel/dfock.h>
#include <src/zcasscf/zsuperci.h>
#include <src/rel/reloverlap.h>

using namespace std;
using namespace bagel;

void ZSuperCI::compute() {
    // TODO : DIIS needs to be introduced eventually

  // ============================
  // macro iteration from here
  // ============================
   Timer timer;

  // intialize coefficients
  init_kramers_coeff();

  cout << setprecision(8) << " kramers restricted re part rms = " << coeff_->get_real_part()->rms() << endl;
  cout << setprecision(8) << " kramers restricted im part rms = " << coeff_->get_imag_part()->rms() << endl;

  if (nact_)
    fci_->update(coeff_);

  for (int iter = 0; iter != 1; ++iter) {

    // first perform CASCI to obtain RDMs
    if (nact_) {
      mute_stdcout(/*fci*/true);
      if (iter) fci_->update(coeff_);
      cout << " Executing FCI calculation in Cycle " << iter << endl;
      fci_->compute();
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      resume_stdcout();
    }


  }
}



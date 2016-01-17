//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcashybrid.cc
// Copyright (C) 2014 Jefferson Bates
//
// Author: Jefferson Bates <jefferson.bates@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <src/multi/zcasscf/zcasscf.h>
#include <src/multi/zcasscf/zsuperci.h>
#include <src/multi/zcasscf/zcasbfgs.h>
#include <src/multi/zcasscf/zcashybrid.h>

 using namespace std;
 using namespace bagel;

void ZCASHybrid::compute() {

  // set convergence threshold
  double global_thresh = idata_->get<double>("thresh", 1.0e-7);
  shared_ptr<Method> active_method;
  // construct and compute SuperCI
  {
    auto idata = make_shared<PTree>(*idata_);
    if (maxiter_switch_ != -1) {
      idata->erase("maxiter");
      idata->put("maxiter", maxiter_switch_);
      idata->put("maxiter_fci", max(maxiter_switch_,50));
    }
    if (thresh_switch_ > 0.0) {
      idata->erase("thresh");
      idata->put("thresh",  thresh_switch_);
    }
    active_method = make_shared<ZSuperCI>(idata, geom_, ref_);
    active_method->compute();
    refout_ = active_method->conv_to_ref();
    complex<double> grad = dynamic_pointer_cast<ZCASSCF>(active_method)->rms_grad();
    if (real(grad) < global_thresh) {
      cout << "      * ZCASSCF converged *    " << endl;
      return;
    }
  }

  // construct and compute step-restricted BFGS
  {
    auto idata = make_shared<PTree>(*idata_);
    idata->erase("active"); // coefficient should be in proper order so active is removed to prevent a second reordering
    idata->erase("kramers_coeff");
    idata->put("kramers_coeff", true); // input coefficient should be time-reversal symmetric by construction
    idata->erase("generate_mvo");
    active_method = make_shared<ZCASBFGS>(idata, geom_, refout_);
    active_method->compute();
    refout_ = active_method->conv_to_ref();
    complex<double> grad = dynamic_pointer_cast<ZCASSCF>(active_method)->rms_grad();
    if (real(grad) < global_thresh) {
      cout << " " << endl;
      cout << "    * ZCASSCF converged *    " << endl;
    } else {
      cout << " " << endl;
      cout << "    * ZCASSCF did NOT converge *    " << endl;
    }
  }

}


shared_ptr<const Reference> ZCASHybrid::conv_to_ref() const {
  assert(refout_);
  return refout_;
}


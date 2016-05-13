//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: casbfgs.h
// Copyright (C) 2014 Toru Shiozaki
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


#include <src/multi/casscf/casscf.h>
#include <src/multi/casscf/superci.h>
#include <src/multi/casscf/casbfgs.h>
#include <src/multi/casscf/cashybrid.h>

 using namespace std;
 using namespace bagel;

void CASHybrid::compute() {

  thresh_ = idata_->get<double>("thresh", 1.0e-8);
  // construct and compute SuperCI
  if (maxiter_switch_ != 0) {
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
    auto active_method = make_shared<SuperCI>(idata, geom_, ref_);
    active_method->compute();
    refout_ = active_method->conv_to_ref();

    const double grad = active_method->rms_grad();
    if (grad < thresh_) {
      cout << "      * CASSCF converged *    " << endl;
      fci_ = active_method->fci();
      energy_ = active_method->energy();
      return;
    }
  }

  // construct and compute step-restricted BFGS
  {
    auto idata = make_shared<PTree>(*idata_);
    idata->erase("active");
    idata->erase("restart");
    auto active_method = make_shared<CASBFGS>(idata, geom_, refout_);
    active_method->compute();
    refout_ = active_method->conv_to_ref();

    const double grad = active_method->rms_grad();
    if (grad < thresh_) {
      cout << " " << endl;
      cout << "      * CASSCF converged *    " << endl;
    }
    fci_ = active_method->fci();
    energy_ = active_method->energy();
  }

}


shared_ptr<const Reference> CASHybrid::conv_to_ref() const {
  assert(refout_);
  return refout_;
}


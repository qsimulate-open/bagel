//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: opt.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#include <functional>
#include <typeinfo>
#include <fstream>
#include <string>
#include <algorithm>
#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/util/io/moldenout.h>
#include <src/wfn/construct_method.h>
#include <src/opt/optimize.h>
#include <src/opt/opt.h>
#include <src/util/archive.h>
#include <src/grad/hess.h>

using namespace std;
using namespace bagel;

// TODO  Constrained optimization

Opt::Opt(shared_ptr<const PTree> idat, shared_ptr<const PTree> inp, shared_ptr<const Geometry> geom, shared_ptr<const Reference> ref)
  : idata_(idat), input_(inp), current_(geom), prev_ref_(ref), iter_(0) {

  auto lastmethod = *idat->get_child("method")->rbegin();
  method_ = to_lower(lastmethod->get<string>("title", ""));

  target_state_ = idat->get<int>("target", 0);
  internal_ = idat->get<bool>("internal", true);
  redundant_ = idat->get<bool>("redundant", false);
  maxiter_ = idat->get<int>("maxiter", 100);
  maxziter_ = idat->get<int>("maxziter", 100);
  scratch_ = idat->get<bool>("scratch", false);
  numerical_ = idat->get<bool>("numerical", false);
  numerical_dx_ = idat->get<double>("numerical_dx", 0.001);
  hess_update_ = idat->get<string>("hess_update", "flowchart");
  hess_approx_ = idat->get<bool>("hess_approx", true);

  constrained_ = idat->get<bool>("constrained", false);
  if (constrained_) {
    if (!internal_ || redundant_) throw runtime_error("Constrained optimization currently only for delocalized internals");
    auto constraints = idat->get_child("constraint");
    for (auto& c : *constraints) {
      constraints_.push_back(make_shared<const OptConstraint>(c));
    }
    cout << "# constraints = " << constraints_.size() << endl;
  }

  explicit_bond_ = idat->get<bool>("explicitbond", false);
  if (explicit_bond_) {
    auto explicit_bonds = idat->get_child("explicit");
    for (auto& e : *explicit_bonds) {
      bonds_.push_back(make_shared<const OptExpBonds>(e));
    }
    cout << endl << "  * Added " << bonds_.size() << " bonds between the non-bonded atoms in overall" << endl;
  }

  if (internal_) {
    if (redundant_)
      bmat_red_ = current_->compute_redundant_coordinate();
    else
      bmat_ = current_->compute_internal_coordinate(nullptr, bonds_, constraints_);
  }

  // small molecule (atomno < 4) threshold : (1.0e-5, 4.0e-5, 1.0e-6)  (tight in GAUSSIAN and Q-Chem = normal / 30)
  // large molecule              threshold : (3.0e-4, 1.2e-3, 1.0e-6)  (normal in GAUSSIAN and Q-Chem)
  if (current_->natom() < 4) {
    thresh_grad_ = idat->get<double>("maxgrad", 0.00001);
    thresh_displ_ = idat->get<double>("maxdisp", 0.00004);
    thresh_echange_ = idat->get<double>("maxchange", 0.000001);
  } else {
    thresh_grad_ = idat->get<double>("maxgrad", 0.0003);
    thresh_displ_ = idat->get<double>("maxdisp", 0.0012);
    thresh_echange_ = idat->get<double>("maxchange", 0.000001);
  }
  opttype_ = idat->get<string>("opttype", "energy");
  maxstep_ = idat->get<double>("maxstep", opttype_ == "energy" ? 0.3 : 0.1);
  algorithm_ = idat->get<string>("algorithm", "ef");
  adaptive_ = idat->get<bool>("adaptive", algorithm_ == "rfo" ? true : false);
#ifndef DISABLE_SERIALIZATION
  refsave_ = idat->get<bool>("save_ref", false);
  if (refsave_) refname_ = idat->get<string>("ref_out", "reference");
#endif

  if (opttype_ == "conical") {
    // parameters for CI optimizations (Bearpark, Robb, Schlegel)
    target_state2_ = idat->get<int>("target2", 1);
    if (target_state2_ > target_state_) {
      int tmpstate = target_state_;
      target_state_ = target_state2_;
      target_state2_ = tmpstate;
    }
    nacmtype_ = idat->get<int>("nacmtype", 1);
    thielc3_  = idat->get<double>("thielc3", 2.0);
    thielc4_  = idat->get<double>("thielc4", 0.5);
    adaptive_ = false;        // we cannot use it for conical intersection optimization because we do not have a target function
  } else if (opttype_ == "mep") {
    // parameters for MEP calculations (Gonzalez, Schlegel)
    mass_weight_ = idat->get<bool>("mass_weight", false);
    mep_direction_ = idat->get<int>("mep_direction", 1);
    if (hess_approx_) throw runtime_error("MEP calculation should be started with Hessian eigenvectors");
  } else if (opttype_ != "energy" && opttype_ != "transition") {
    throw runtime_error("Optimization type should be: \"energy\", \"transition\", \"conical\", or \"mep\"");
  }

}

void Opt::compute() {
  auto displ = make_shared<XYZFile>(current_->natom());
  size_ = internal_ ? (redundant_? bmat_red_[0]->ndim() : bmat_[0]->mdim()) : current_->natom()*3;

  dispsize_ = max(current_->natom(),int(size_/3+1));
  displ_ = make_shared<XYZFile>(dispsize_);
  grad_ = make_shared<GradFile>(dispsize_);

  auto mep_start = make_shared<XYZFile>(current_->natom());

  if (hess_approx_) {
    cout << "    * Use approximate Hessian for optimization" << endl;
    if (internal_ && !redundant_) hess_ = make_shared<Matrix>(*(bmat_[2]));
    else {
      hess_ = make_shared<Matrix>(size_, size_);
      hess_->unit();
    }
  } else {
    cout << "    * Compute molecular Hessian for optimization" << endl;
    auto hess = make_shared<Hess>(idata_, current_, prev_ref_);
    hess->compute();
    // if internal, we should transform the Hessian according to Eq. (6) in Schlegel
    // dB/dX term currently omitted (reasonable approximation: Handbook of Computational Chemistry, pp. 323--324)
    if (internal_)
      if (redundant_)
        hess_ = make_shared<Matrix>(*bmat_red_[1] % *(hess->hess()) * *bmat_red_[1]);
      else
        hess_ = make_shared<Matrix>(*bmat_[1] % *(hess->hess()) * *bmat_[1]);
    else
      hess_ = hess->hess()->copy();

    copy_n(hess->proj_hess()->element_ptr(0,abs(mep_direction_) - 1), current_->natom()*3, mep_start->data());
  }

  if (opttype_ == "mep") {
    do_mep(mep_start);
  } else {
    do_optimize();
  }

#ifndef DISABLE_SERIALIZATION
  if (refsave_) {
    OArchive archive(refname_);
    archive << prev_ref_;
  }
#endif
}


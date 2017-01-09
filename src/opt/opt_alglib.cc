//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: opt_alglib.cc
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

using namespace std;
using namespace bagel;

void Opt::compute_alglib() {
  assert(typeid(double) == typeid(double));
  std::shared_ptr<const XYZFile> displ = current_->xyz();
  size_ = internal_ ? bmat_[0]->mdim() : displ->size();

  if (internal_)
    displ = displ->transform(bmat_[0], true);

  try {
    alglib::real_1d_array x;
    x.setcontent(size_, displ->data());
    eval_type eval = std::bind(&Opt::evaluate_alglib, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
  
    if (algorithm_ == "cg") {
      alglib::mincgstate state;
      alglib::mincgreport rep;
  
      alglib::mincgcreate(x, state);
      alglib::mincgsetcond(state, thresh_*std::sqrt(size_), 0.0, 0.0, maxiter_);
      alglib::mincgsetstpmax(state, maxstep_);
  
      alglib::mincgoptimize(state, eval);
  
    } else if (algorithm_ == "lbfgs") {
      alglib::minlbfgsstate state;
      alglib::minlbfgsreport rep;
  
      alglib::minlbfgscreate(1, x, state); 
      alglib::minlbfgssetcond(state, thresh_*std::sqrt(size_), 0.0, 0.0, maxiter_);
      alglib::minlbfgssetstpmax(state, maxstep_);
  
      alglib::minlbfgsoptimize(state, eval);
  
    }
  } catch (alglib::ap_error e) {
    std::cout << e.msg << std::endl;
    throw std::runtime_error("optimization failed");
  }
}

void Opt::evaluate_alglib(const alglib::real_1d_array& x, double& en, alglib::real_1d_array& grad, void* ptr) {
  std::shared_ptr<const XYZFile> xyz = current_->xyz();

  // first convert x to the geometry
  auto displ = std::make_shared<XYZFile>(current_->natom());
  assert(size_ == x.length());
  std::copy_n(x.getcontent(), size_, displ->data()); 

  if (internal_)
    displ = displ->transform(bmat_[1], false);

  *displ -= *xyz;

  muffle_ = make_shared<Muffle>("opt.log");
  if (iter_ > 0) muffle_->mute();

  // current Geometry
  if (iter_ > 0) {
    current_ = std::make_shared<Geometry>(*current_, displ, std::make_shared<const PTree>()); 
    current_->print_atoms();
    if (internal_)
      bmat_ = current_->compute_internal_coordinate(bmat_[0]);
  }

  // first calculate reference (if needed)
  std::shared_ptr<PTree> cinput; 
  std::shared_ptr<const Reference> ref;
  if (!prev_ref_ || scratch_) {
    auto m = input_->begin();
    for ( ; m != --input_->end(); ++m) {
      const std::string title = to_lower((*m)->get<std::string>("title", ""));
      if (title != "molecule") {
        std::shared_ptr<Method> c = construct_method(title, *m, current_, ref);
        if (!c) throw std::runtime_error("unknown method in optimization");
        c->compute();
        ref = c->conv_to_ref();
      } else {
        current_ = std::make_shared<const Geometry>(*current_, *m); 
        if (ref) ref = ref->project_coeff(current_);
      }
    }
    cinput = std::make_shared<PTree>(**m);
  } else {
    ref = prev_ref_->project_coeff(current_);
    cinput = std::make_shared<PTree>(**input_->rbegin());
  }
  cinput->put("gradient", true);

  // then calculate gradients
  double rms;
  {
    // current geom and grad in the cartesian coordinate
    if (iter_ == 0) {
      print_header();
      muffle_->mute();
    }
    shared_ptr<GradFile> cgrad = get_grad(cinput, ref);
    en = en_;
    if (internal_)
      cgrad = cgrad->transform(bmat_[1], true);

    assert(size_ == grad.length());
    std::copy_n(cgrad->data(), size_, grad.getcontent());

    // current geometry in a molden file
    MoldenOut mfs("opt.molden");
    mfs << current_;
    rms = cgrad->rms();
  }

  ++iter_;
  // returns energy

  muffle_->unmute();
  print_iteration(rms, timer_.tick());
}

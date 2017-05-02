//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: do_opt.cc
// Copyright (C) 2017 Toru Shiozaki
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

void Opt::do_optimize() {
  auto displ = make_shared<XYZFile>(current_->natom());

  muffle_ = make_shared<Muffle>("opt.log");
  muffle_->unmute();

  for (iter_ = 1; iter_ != maxiter_; ++iter_) {
    shared_ptr<const XYZFile> xyz = current_->xyz(); 
    prev_xyz_.push_back(xyz);

    displ = displ_;

    if (internal_) {
      if (redundant_)
        displ = displ->transform(bmat_red_[1], false);
      else
        displ = iterate_displ();
    }

    prev_displ_.push_back(displ);
    current_ = make_shared<Geometry>(*current_, displ, make_shared<const PTree>());
    current_->print_atoms();
    if (internal_) {
      if (redundant_)
        bmat_red_ = current_->compute_redundant_coordinate(bmat_red_[0]);
      else
        bmat_ = current_->compute_internal_coordinate(bmat_[0], bonds_, constraints_);
    }

    shared_ptr<PTree> cinput;
    shared_ptr<const Reference> ref;
    if (!prev_ref_ || scratch_) {
      auto m = input_->begin();
      for ( ; m != --input_->end(); ++m) {
        const string title = to_lower((*m)->get<string>("title", ""));
        if (title != "molecule") {
          shared_ptr<Method> c = construct_method(title, *m, current_, ref);
          if (!c) throw runtime_error("unknown method in optimization");
          c->compute();
          ref = c->conv_to_ref();
        } else {
          current_ = make_shared<const Geometry>(*current_, *m);
          if (ref) ref = ref->project_coeff(current_);
        }
      }
      cinput = make_shared<PTree>(**m);
    } else {
      ref = prev_ref_->project_coeff(current_);
      cinput = make_shared<PTree>(**input_->rbegin());
    }
    cinput->put("gradient_", true);

    double rms;
    double maxgrad;
    {
      grad_->zero();
      shared_ptr<GradFile> cgrad = get_grad(cinput, ref);
      prev_grad_.push_back(cgrad);
      grad_->add_block(1.0, 0, 0, 3, current_->natom(), cgrad);
      rms = cgrad->rms();       // This is more appropriate
      maxgrad = cgrad->maximum(current_->natom());

      if (internal_) {
        if (redundant_)
          grad_ = grad_->transform(bmat_red_[1], true);
        else
          grad_ = grad_->transform(bmat_[1], true);
      }

      // Update Hessian with Flowchart method
      if (iter_ != 1)
        hessian_update();

      prev_grad_internal_ = make_shared<GradFile>(*grad_);

      MoldenOut mfs("opt.molden");
      mfs << current_;

      displ_ = get_step();
    }

    displ = displ_;

    // check the size of (displ)
    if (internal_) {
      if (redundant_)
        displ = displ->transform(bmat_red_[1], false);
      else
        displ = displ->transform(bmat_[1], false);
    }

    if (adaptive_) do_adaptive();
    double maxdispl = displ->maximum(current_->natom());
    double echange = en_ - (prev_en_.empty() ? 0.0 : prev_en_.back());

    bool convergegrad = maxgrad > thresh_grad_ ? false : true;
    bool convergedispl = maxdispl > thresh_displ_ ? false : true;
    bool convergeenergy = fabs(echange) > thresh_echange_ ? false : true;
    cout << endl << "  === Convergence status ===" << endl << endl;
    cout << "                         Maximum     Tolerance   Converged?" << endl;
    cout << "  * Gradient      " << setw(14) << setprecision(6) << maxgrad << setw(14) << thresh_grad_ << setw(13) << (convergegrad? "Yes" : "No") << endl;
    cout << "  * Displacement  " << setw(14) << setprecision(6) << maxdispl << setw(14) << thresh_displ_ << setw(13) << (convergedispl? "Yes" : "No") << endl;
    cout << "  * Energy change " << setw(14) << setprecision(6) << echange << setw(14) << thresh_echange_ << setw(13) << (convergeenergy? "Yes" : "No") << endl << endl;

    muffle_->unmute();
    if (iter_ == 1)
      print_header();
    prev_en_.push_back(en_);
    print_iteration(rms, timer_.tick());

    muffle_->mute();

    if (convergegrad && (convergedispl || convergeenergy)) break;
  }
}

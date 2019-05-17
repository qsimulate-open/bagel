//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: compute_opt.cc
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
#include <src/grad/hess.h>
#include <src/grad/gradeval.h>
#include <src/util/timer.h>
#include <src/util/io/moldenout.h>
#include <src/opt/optimize.h>
#include <src/opt/opt.h>

using namespace std;
using namespace bagel;

void Opt::compute_optimize() {
  auto displ = make_shared<XYZFile>(current_->natom());

  muffle_ = make_shared<Muffle>("opt.log");
  muffle_->unmute();

  for (iter_ = 1; iter_ != optinfo()->maxiter(); ++iter_) {
    const bool hess_recalc = (iter_ % optinfo()->hess_recalc_freq() == 1 && iter_ != 1 && optinfo()->hess_approx() == false && optinfo()->opttype()->is_transition());
    shared_ptr<const XYZFile> xyz = current_->xyz();
    prev_xyz_.push_back(xyz);

    displ = displ_;

    if (optinfo()->internal() && iter_ != 1) {
      if (optinfo()->redundant()) {
        displ = displ->transform(bmat_red_[1], false);
      } else {
        displ = iterate_displ();
      }
    }

    if (iter_ != 1) {
      prev_displ_.push_back(displ);
      current_ = make_shared<Geometry>(*current_, displ, make_shared<const PTree>());
      current_->print_atoms();
    }

    if (optinfo()->internal()) {
      if (optinfo()->redundant())
        tie(bondlist_, bmat_red_) = current_->compute_redundant_coordinate(bondlist_);
      else
        bmat_ = current_->compute_internal_coordinate(bmat_[0], optinfo()->bonds(), optinfo()->opttype()->is_transition(), false);
    }

    shared_ptr<PTree> cinput;
    shared_ptr<const Reference> ref;
    tie(cinput, ref, current_) = get_grad_input();

    double rms;
    double maxgrad;
    double param;
    {
      grad_->zero();

      shared_ptr<GradFile> cgrad;
      tie(en_, param, prev_ref_, cgrad) = get_grad(cinput, ref);
      prev_grad_.push_back(cgrad);
      grad_->add_block(1.0, 0, 0, 3, current_->natom(), cgrad);

      rms = cgrad->rms();
      maxgrad = cgrad->maximum(current_->natom());

      if (optinfo()->internal()) {
        if (optinfo()->redundant()) {
          grad_ = grad_->transform(bmat_red_[1], true);
          prev_xyz_internal_.push_back(make_shared<GradFile>(*(current_->xyz()->transform(bmat_red_[0], true))));
        } else {
          grad_ = grad_->transform(bmat_[1], true);
          prev_xyz_internal_.push_back(make_shared<GradFile>(*(current_->xyz()->transform(bmat_[0], true))));
        }
      } else {
        prev_xyz_internal_.push_back(make_shared<XYZFile>(*xyz));
      }
      prev_grad_internal_.push_back(make_shared<GradFile>(*grad_));


      // Update Hessian with Flowchart method
      if (hess_recalc) {
        auto hess = make_shared<Hess>(idata_, current_, prev_ref_);
        hess->compute();
        if (optinfo()->internal()) {
          if (optinfo()->redundant())
            hess_ = make_shared<Matrix>(*bmat_red_[1] % *(hess->hess()) * *bmat_red_[1]);
          else
            hess_ = make_shared<Matrix>(*bmat_[1] % *(hess->hess()) * *bmat_[1]);
        } else {
          hess_ = hess->hess()->copy();
        }
      } else if (iter_ != 1 && !(optinfo()->hess_approx() == false && optinfo()->opttype()->is_transition() && iter_ % optinfo()->hess_recalc_freq() == 2)) {
        hess_ = hessian_update();
      }

      MoldenOut mfs("opt.molden");
      mfs << current_;

      tie(predictedchange_, predictedchange_prev_, displ_) = get_step();

    }

    displ = displ_;

    // check the size of (displ)
    if (optinfo()->internal()) {
      if (optinfo()->redundant())
        displ = displ->transform(bmat_red_[1], false);
      else
        displ = displ->transform(bmat_[1], false);
    }

    if (optinfo()->adaptive())
      maxstep_ = do_adaptive(iter_);

    const double maxdispl = displ->maximum(current_->natom());
    const double echange = en_ - (prev_en_.empty() ? 0.0 : prev_en_.back());

    const bool convergegrad = !(maxgrad > optinfo()->thresh_grad());
    const bool convergedispl = !(maxdispl > optinfo()->thresh_displ());
    const bool convergeenergy = !(fabs(echange) > optinfo()->thresh_echange());

    cout << endl << "  === Convergence status ===" << endl << endl;
    cout << "                         Maximum     Tolerance   Converged?" << endl;
    cout << "  * Gradient      " << setw(14) << setprecision(6) << maxgrad << setw(14) << optinfo()->thresh_grad() << setw(13) << (convergegrad? "Yes" : "No") << endl;
    cout << "  * Displacement  " << setw(14) << setprecision(6) << maxdispl << setw(14) << optinfo()->thresh_displ() << setw(13) << (convergedispl? "Yes" : "No") << endl;
    cout << "  * Energy change " << setw(14) << setprecision(6) << echange << setw(14) << optinfo()->thresh_echange() << setw(13) << (convergeenergy? "Yes" : "No") << endl << endl;

    muffle_->unmute();

    if (iter_ == 1)
      print_header();

    prev_en_.push_back(en_);
    print_iteration(iter_, rms, param, timer_.tick());

    muffle_->mute();

    if (convergegrad && (convergedispl || convergeenergy))
      break;
  }

  muffle_->unmute();
}

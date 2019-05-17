//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: compute_mep.cc
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
#include <src/opt/optimize.h>
#include <src/opt/opt.h>
#include <src/util/archive.h>

using namespace std;
using namespace bagel;

void Opt::compute_mep(shared_ptr<XYZFile> mep_start) {
  // performs second order MEP calculation in Cartesian or in internal coordinates (J. Chem. Phys. 1989, 90, 2154)
  cout << "    * Doing second order MEP calculation" << endl;

  if (optinfo()->mep_direction() < 0)
    mep_start->scale(-1.0);
  else
    mep_start->scale(1.0);

  if (optinfo()->internal()) {
    if (optinfo()->redundant())
      mep_start = mep_start->transform(bmat_red_[1], true);
    else
      mep_start = mep_start->transform(bmat_[1], true);
  }

  muffle_ = make_shared<Muffle>("mep.log");
  muffle_->mute();

  // get gradient at the current point
  {
    shared_ptr<PTree> cinput;
    shared_ptr<const Reference> ref;
    tie(cinput, ref, current_) = get_grad_input();

    {
      grad_->zero();
      shared_ptr<GradFile> cgrad;
      tie(en_, ignore, prev_ref_, cgrad) = get_grad(cinput, ref);
      grad_->add_block(1.0, 0, 0, 3, current_->natom(), cgrad);

      if (optinfo()->internal()) {
        if (optinfo()->redundant())
          grad_ = grad_->transform(bmat_red_[1], true);
        else
          grad_ = grad_->transform(bmat_[1], true);
      }

      if (optinfo()->mep_direction() == 0) {
        copy_n(grad_->data(), size_, mep_start->data());
        mep_start->scale(-1.0);
      }
    }
  }

  // get step
  mep_start->scale(maxstep_ / mep_start->norm());
  copy_n(mep_start->data(), size_, displ_->data());

  muffle_->unmute();

  // macroiteration
  auto displ = make_shared<XYZFile>(current_->natom());

  for (int iter = 1; iter != optinfo()->maxiter(); ++iter) {
    if (optinfo()->internal()) {
      if (optinfo()->redundant()) {
        grad_ = grad_->transform(bmat_red_[1], true);
        prev_xyz_internal_.push_back(make_shared<GradFile>(*(current_->xyz()->transform(bmat_red_[0], true))));
      } else {
        grad_ = grad_->transform(bmat_[1], true);
        prev_xyz_internal_.push_back(make_shared<GradFile>(*(current_->xyz()->transform(bmat_[0], true))));
      }
    } else {
      prev_xyz_internal_.push_back(make_shared<XYZFile>(*(current_->xyz())));
    }
    prev_grad_internal_.push_back(make_shared<GradFile>(*grad_));

    prev_en_.push_back(en_);
    prev_xyz_.push_back(current_->xyz());
    cout << endl << " ============================ MEP point # " << setw(4) << iter << " ============================" << endl;
    current_->print_atoms();
    cout << "    MEP energy at # " << setw(4) << iter << " : " << setprecision(10) << en_ << endl;
    cout << " ===========================================================================" << endl << endl;
    muffle_->mute();
    if (iter != 1) {
      copy_n(grad_->data(), size_, mep_start->data());
      mep_start->scale(maxstep_ / mep_start->norm());
      copy_n(mep_start->data(), size_, displ_->data());
      displ_->scale(-1.0);
    }

    // x_{k+1}^{\*} = x_k + 0.5 * p
    displ_->scale(0.5);
    displ = displ_;
    if (optinfo()->internal() && iter != 1) {
      if (optinfo()->redundant())
        displ = displ->transform(bmat_red_[1], false);
      else
        displ = iterate_displ();
    }

    if (iter != 1) {
      current_ = make_shared<Geometry>(*current_, displ, make_shared<const PTree>());
      current_->print_atoms();
    }

    if (optinfo()->internal()) {
      if (optinfo()->redundant())
        tie(bondlist_, bmat_red_) = current_->compute_redundant_coordinate(bondlist_);
      else
        bmat_ = current_->compute_internal_coordinate(bmat_[0], optinfo()->bonds(), false, false);
    }

    if (iter > 1) {
      // Hessian updater
      auto y  = make_shared<GradFile>(*(prev_grad_internal_[iter-1]) - *(prev_grad_internal_[iter-2]));
      auto s  = make_shared<GradFile>((*prev_xyz_internal_[iter-1] - *prev_xyz_internal_[iter-2]));
      auto hs = make_shared<GradFile>(*(s->transform(hess_, /*transpose=*/false)));
//      auto dhess = hessian_update_bfgs(y, s, hs);
      auto dhess = hessian_update_bofill(y, s, hs);
      *hess_ += *dhess;
    }

    // microiteration with quasi-Newton search
    // initial guess: p_{k+1} = x_{k+1}^{\prime} - x_k
    auto pd = make_shared<XYZFile>(current_->natom());
    pd = displ_;
    bool flag = false;
    for (int miciter = 1; miciter != optinfo()->maxiter(); ++miciter) {

      // move geometry
      auto dx = make_shared<XYZFile>(current_->natom());
      dx = displ_;
      if (optinfo()->internal()) {
        if (optinfo()->redundant())
          dx = dx->transform(bmat_red_[1], false);
        else
          dx = iterate_displ();
      }

      current_ = make_shared<Geometry>(*current_, dx, make_shared<const PTree>());
      current_->print_atoms();
      if (optinfo()->internal()) {
        if (optinfo()->redundant())
          tie(bondlist_, bmat_red_) = current_->compute_redundant_coordinate(bondlist_);
        else
          bmat_ = current_->compute_internal_coordinate(bmat_[0], optinfo()->bonds(), false, false);
      }

      // compute gradient
      shared_ptr<PTree> cinput;
      shared_ptr<const Reference> ref;
      tie(cinput, ref, current_) = get_grad_input();

      {
        grad_->zero();
        shared_ptr<GradFile> cgrad;
        tie(en_, ignore, prev_ref_, cgrad) = get_grad(cinput, ref);
        grad_->add_block(1.0, 0, 0, 3, current_->natom(), cgrad);

        if (optinfo()->internal()) {
          if (optinfo()->redundant())
            grad_ = grad_->transform(bmat_red_[1], true);
          else
            grad_ = grad_->transform(bmat_[1], true);
        }
      }

      {
        auto vj = make_shared<XYZFile>(dispsize_);
        VectorB pj(size_);
        VectorB gj(size_);
        auto hess = make_shared<Matrix>(*hess_);
        VectorB bj(size_);
        hess->diagonalize(bj);

        for (int i = 0; i != size_; ++i) {
          copy_n(hess->element_ptr(0,i), size_, vj->data());
          pj[i] = vj->dot_product(pd);
          gj[i] = vj->dot_product(grad_);
        }

        // Newton-Raphson search for lambda
        // initial value: lambda should be smaller than the smallest eigenvalue
        cout << endl << endl << "  === lambda iteration ===" << endl;
        double lambda = bj[0] - 0.1;
        for (int iiter = 0; iiter != 100; ++iiter) {
          const double lambda_prev = lambda;
          double fv = 0.0;
          double df = 0.0;
          for (int i = 0; i != size_; ++i) {
            const double bpgb = (bj[i] * pj[i] - gj[i]) / (bj[i] - lambda);
            const double bl   = bj[i] - lambda;
            fv += bpgb * bpgb;
            df += 2.0 * bpgb * bpgb / bl;
          }
          fv -= 0.25 * maxstep_ * maxstep_;
          lambda -= fv / df;
          cout << setw(7) << iiter << setw(15) << setprecision(10) << fabs(fv) << endl;
          if (fabs(lambda - lambda_prev) < 1.0e-8)
            break;
          if (iiter == 99)
            flag = true;
        }
        cout << endl << endl;

        dx = make_shared<XYZFile>(*grad_ - *pd * lambda);
        auto hinv = make_shared<Matrix>(*hess_);
        auto unit = make_shared<Matrix>(*hess_);
        unit->unit();
        *hinv = *hinv - *unit * lambda;
        hinv->inverse();
        dx = dx->transform(hinv, /*transpose=*/false);
        dx->scale(-1.0);
      }
      if (flag)
        break;

      *pd = *pd + *dx;
      displ_ = dx;

      auto tangent_vec = make_shared<GradFile>(*grad_ - *pd * (pd->dot_product(*grad_) / (pd->norm() * pd->norm())));
      cout << "  * MEP microiteration " << miciter << " residual = " << setprecision(10) << tangent_vec->norm() << "  " << dx->norm() << endl;
      if (tangent_vec->norm() < 1.0e-5 && dx->norm() < 1.0e-5)
        break;

    }
    muffle_->unmute();

    if (fabs(prev_en_.back() - en_) < 1.0e-6 && !(optinfo()->mep_direction() == 0 && iter == 1)) {
      cout << "  * MEP job converged to reactant / product, or too small stepsize." << endl;
      break;
    } else if (flag) {
      cout << "  * MEP job encountered convergence problem in lambda iteration." << endl << "    Possibly converged near reactant / product." << endl;
      break;
    }
  }
}

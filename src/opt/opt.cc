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

using namespace std;
using namespace bagel;

// TODO: (1) Introduce redundant internal coordinates
//       (2) Add eigenvector following algorithm
//       (3) Scaled RFO algorithm (different from rfos)
//       (4) Transition state search

void Opt::compute_noalglib() {
  auto displ = make_shared<XYZFile>(current_->natom());
  size_ = internal_ ? (redundant_? bmat_red_[0]->ndim() : bmat_[0]->mdim()) : displ->size();
  
  dispsize_ = max(current_->natom(),int(size_/3+1));
  displ_ = make_shared<XYZFile>(dispsize_);
  grad_ = make_shared<GradFile>(dispsize_);

  hess_ = make_shared<Matrix>(size_, size_);
  hess_->unit();        // TODO can take the initial hessian from internal coordinate generator
  
  print_header();

  for (iter_ = 1; iter_ != maxiter_; ++iter_) {
    shared_ptr<const XYZFile> xyz = current_->xyz();
    mute_stdcout();

    displ = displ_;

    if (internal_) {
      if (redundant_) 
        displ = displ->transform(bmat_red_[2], false);      // should be done iteratively
      else 
        displ = displ->transform(bmat_[1], false);
    }

    current_ = make_shared<Geometry>(*current_, displ, make_shared<const PTree>()); 
    current_->print_atoms();
    if (internal_) {
      if (redundant_) 
        bmat_red_ = current_->compute_redundant_coordinate(bmat_red_[0]);
      else            
        bmat_ = current_->compute_internal_coordinate(bmat_[0]);
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
    cinput->put("gradient", true);
 
    double rms;
    {
      grad_->zero();
      shared_ptr<GradFile> cgrad = get_grad(cinput, ref);
      grad_->add_block(1.0, 0, 0, 3, current_->natom(), cgrad);

      if (internal_) {
        if (redundant_)
          grad_ = grad_->transform(bmat_red_[1], false);
        else 
          grad_ = grad_->transform(bmat_[1], true);
      }

      // Update Hessian with Flowchart method
      if (iter_ != 1)
        hessian_update();

      prev_grad_ = make_shared<GradFile>(*grad_);
 
      MoldenOut mfs("opt.molden");
      mfs << current_;
      rms = grad_->rms();

      displ_ = get_step();
    }
    
    if (adaptive_) do_adaptive();

    resume_stdcout();
    print_iteration(rms, timer_.tick());
    en_prev_ = en_;

    double stepnorm = displ_->norm();
    if (stepnorm < thresh_ * sqrt(size_)) break;
  }

#ifndef DISABLE_SERIALIZATION
  if (refsave_) {
    OArchive archive(refname_);
    archive << prev_ref_;
  }
#endif
}

void Opt::hessian_update() {
  auto y  = make_shared<GradFile>(*grad_ - *prev_grad_);
  auto s  = make_shared<GradFile>(*displ_);
  auto hs = make_shared<GradFile>(*(s->transform(hess_, /*transpose=*/false)));
  auto z  = make_shared<GradFile>(*y - *hs);

  double nzs = z->norm() * s->norm();
  double nys = y->norm() * s->norm();
  double zs = z->dot_product(s);
  double ys = y->dot_product(s);

  if ((zs / nzs) < -0.1) hessian_update_sr1(y,s,z);
  else if ((ys / nys) > 0.1) hessian_update_bfgs(y,s,hs);
  else hessian_update_psb(y,s,z);
}

void Opt::hessian_update_sr1(shared_ptr<GradFile> y, shared_ptr<GradFile> s, shared_ptr<GradFile> z) {
  cout << "Updating Hessian using SR1 " << endl;
  // Hessian update with SR1
  double  zs = z->dot_product(s);
  if (fabs(zs)>1.0e-12) zs = 1.0 / zs;

  auto zzt = make_shared<Matrix>(size_,size_);
  dger_(size_,size_,zs,z->data(),1,z->data(),1,zzt->data(),size_);

  *hess_ = *hess_ + *zzt;
}

void Opt::hessian_update_bfgs(shared_ptr<GradFile> y, shared_ptr<GradFile> s, shared_ptr<GradFile> hs) {
  cout << "Updating Hessian using BFGS " << endl;
  // Hessian update with BFGS
  double shs = hs->dot_product(s);
  double  ys = y->dot_product(s);
  if (fabs(ys)>1.0e-12) ys = 1.0 / ys;
  if (fabs(shs)>1.0e-12) shs = -1.0 / shs;

  auto yyt = make_shared<Matrix>(size_,size_);
  auto sst = make_shared<Matrix>(size_,size_);
  dger_(size_,size_,shs,s->data(),1,s->data(),1,sst->data(),size_);
  dger_(size_,size_,ys,y->data(),1,y->data(),1,yyt->data(),size_);

  auto bsst = make_shared<Matrix>(*hess_ * *sst * *hess_);

  *hess_ = *hess_ + *bsst + *yyt;
}

void Opt::hessian_update_psb(shared_ptr<GradFile> y, shared_ptr<GradFile> s, shared_ptr<GradFile> z) {
  cout << "Updating Hessian using PSB " << endl;
  // Hessian update with PSB
  double ss = s->dot_product(s);
  double ss2 = ss * ss;
  double sz = s->dot_product(z);

  ss = 1.0 / ss;
  ss2= -sz / ss2;

  auto szt = make_shared<Matrix>(size_,size_);
  auto zst = make_shared<Matrix>(size_,size_);
  auto sst = make_shared<Matrix>(size_,size_);
  dger_(size_,size_,ss,s->data(),1,z->data(),1,szt->data(),size_);
  dger_(size_,size_,ss,z->data(),1,s->data(),1,zst->data(),size_);
  dger_(size_,size_,ss2,s->data(),1,s->data(),1,sst->data(),size_);

  *hess_ = *hess_ + *szt + *zst + *sst;
}

shared_ptr<XYZFile> Opt::get_step() {
  auto displ = make_shared<XYZFile>(dispsize_);

  if (algorithm_ == "steep")
    displ = get_step_steep();
  else if (algorithm_ == "nr")
    displ = get_step_nr();
  else if (algorithm_ == "rfo")
    displ = get_step_rfo();
  else if (algorithm_ == "rfos")
    displ = get_step_rfos();

  return displ;
}

shared_ptr<XYZFile> Opt::get_step_steep() {
  // Steepest descent step
  auto displ = make_shared<XYZFile>(dispsize_);
  copy_n(grad_->data(), size_, displ->data());
  displ->scale(-1.0);

  double stepnorm = displ->norm();
  if (stepnorm > (maxstep_))
    displ->scale(maxstep_/stepnorm);

  return displ;
}

shared_ptr<XYZFile> Opt::get_step_nr() {

  // Quasi-Newton-Raphson step
  auto displ = make_shared<XYZFile>(dispsize_);

  shared_ptr<Matrix> hinv(hess_);
  hinv->inverse();

  copy_n(grad_->data(), size_, displ->data());
  displ = displ->transform(hinv, /*transpose=*/false);

  displ->scale(-1.0);

  double stepnorm = displ->norm();
  if (stepnorm > (maxstep_))
    displ->scale(maxstep_/stepnorm);

  return displ;
}

shared_ptr<XYZFile> Opt::get_step_rfo() {
  // Rational function optimization (aka augmented Hessian)
  // Here we do scale lambda to get step < steplength

  auto displ = make_shared<XYZFile>(dispsize_);
  {
    auto aughes = make_shared<Matrix>(size_+1,size_+1);
    VectorB eigv(size_+1);
    double lambda = 1.0;
    while (1) {
      aughes->zero();
      aughes->add_block(1.0 * lambda, 1, 1, size_, size_, hess_);
      aughes->add_block(1.0, 1, 0, size_, 1, grad_->data());
      aughes->add_block(1.0, 0, 1, 1, size_, grad_->data());
  
      aughes->diagonalize(eigv);
      aughes->scale(lambda / aughes->element(0,0));
 
      copy_n(aughes->element_ptr(1,0), size_, displ->data());
      if (displ->norm() < maxstep_) break;
      else lambda /= 1.2;
    }
  }

  if (adaptive_) {
    // When we use adaptive steplength, we should predict quadratic energy change

    double qg  = displ->dot_product(grad_);
    auto hq = make_shared<GradFile>(*(displ->transform(hess_, /*transpose=*/false)));
    double qhq = displ->dot_product(hq);
    predictedchange_prev_ = predictedchange_;
    predictedchange_ = qg + 0.5 * qhq;
  }

  return displ;
}

shared_ptr<XYZFile> Opt::get_step_rfos() {
  // Rational function optimization (aka augmented Hessian)
  // Here we do not scale lambda, but just scale the computed step

  auto displ = make_shared<XYZFile>(dispsize_);
  {
    auto aughes = make_shared<Matrix>(size_+1,size_+1);
    VectorB eigv(size_+1);
    aughes->zero();
    aughes->add_block(1.0, 1, 1, size_, size_, hess_);
    aughes->add_block(1.0, 1, 0, size_, 1, grad_->data());
    aughes->add_block(1.0, 0, 1, 1, size_, grad_->data());
  
    aughes->diagonalize(eigv);
    aughes->scale(1.0 / aughes->element(0,0));
 
    copy_n(aughes->element_ptr(1,0), size_, displ->data());
    if (displ->norm() > maxstep_) displ->scale(maxstep_ / displ->norm());
  }

  if (adaptive_) {
    // When we use adaptive steplength, we should predict quadratic energy change

    double qg  = displ->dot_product(grad_);
    auto hq = make_shared<GradFile>(*(displ->transform(hess_, /*transpose=*/false)));
    double qhq = displ->dot_product(hq);
    predictedchange_prev_ = predictedchange_;
    predictedchange_ = qg + 0.5 * qhq;
  }

  return displ;
}

void Opt::do_adaptive() {
  // Fletcher's adaptive stepsize algorithm

  bool algo = iter_ > 1 && (algorithm_ == "rfo" || algorithm_ == "rfos");

  if (algo) {
    realchange_ = en_ - en_prev_;
    double predreal_ratio = realchange_ / predictedchange_prev_;
    if (predreal_ratio > 1.0) predreal_ratio = 1.0 / predreal_ratio;
    
    if (predreal_ratio > 0.75 && displ_->norm() > 0.80*maxstep_) maxstep_ *= 2.0;
    else if (predreal_ratio < 0.25) maxstep_ *= 0.25;
  }
}

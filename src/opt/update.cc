//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: update.cc
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

shared_ptr<Matrix> Opt::hessian_update() const {
  auto y  = make_shared<GradFile>(*grad_ - *(prev_grad_internal_[iter_-2]));
  auto s  = make_shared<GradFile>(*displ_);
  auto hs = make_shared<GradFile>(*(s->transform(hess_, /*transpose=*/false)));
  auto z  = make_shared<GradFile>(*y - *hs);
  shared_ptr<Matrix> dhess;

  if (optinfo()->hessupdate()->is_flowchart()) {

    const double nzs = z->norm() * s->norm();
    const double nys = y->norm() * s->norm();
    const double zs = z->dot_product(s);
    const double ys = y->dot_product(s);

    if ((zs / nzs) < -0.1) {
      cout << "  * Updating Hessian using SR1 " << endl;
      dhess = hessian_update_sr1(y, s, z);
    } else if ((ys / nys) > 0.1) {
      cout << "  * Updating Hessian using BFGS " << endl;
      dhess = hessian_update_bfgs(y, s, hs);
    } else {
      cout << "  * Updating Hessian using PSB " << endl;
      dhess = hessian_update_psb(y, s, z);
    }

  } else if (optinfo()->hessupdate()->is_bfgs()) {

    cout << "  * Updating Hessian using BFGS " << endl;
    dhess = hessian_update_bfgs(y, s, hs);

  } else if (optinfo()->hessupdate()->is_psb()) {

    cout << "  * Updating Hessian using PSB " << endl;
    dhess = hessian_update_psb(y, s, z);

  } else if (optinfo()->hessupdate()->is_sr1()) {

    cout << "  * Updating Hessian using SR1 " << endl;
    dhess = hessian_update_sr1(y, s, z);

  } else if (optinfo()->hessupdate()->is_bofill()) {

    cout << "  * Updating Hessian using Bofill's algorithm for combining SR1 and PSB" << endl;
    dhess = hessian_update_bofill(y, s, z);

  } else if (optinfo()->hessupdate()->is_noupdate()) {

    cout << "  * Does not update Hessian" << endl;
    dhess = hess_->clone();

  } else {

    throw runtime_error ("available Hessian update schemes are: \"flowchart\", \"bfgs\", \"psb\", \"sr1\", and \"bofill\"");
  }

  *dhess += *hess_;

  return dhess;
}


shared_ptr<Matrix> Opt::hessian_update_sr1(shared_ptr<const GradFile> y, shared_ptr<const GradFile> s, shared_ptr<const GradFile> z) const {
  // Hessian update with SR1

  double  zs = z->dot_product(s);
  if (fabs(zs) > 1.0e-12)
    zs = 1.0 / zs;

  auto zzt = make_shared<Matrix>(size_, size_);
  dger_(size_, size_, zs, z->data(), 1, z->data(), 1, zzt->data(), size_);

  auto hess = make_shared<Matrix>(*zzt);

  return hess;
}


shared_ptr<Matrix> Opt::hessian_update_bfgs(shared_ptr<const GradFile> y, shared_ptr<const GradFile> s, shared_ptr<const GradFile> hs) const {
  // Hessian update with BFGS
  double shs = hs->dot_product(s);
  double  ys = y->dot_product(s);
  if (fabs(ys) > 1.0e-12)
    ys = 1.0 / ys;
  if (fabs(shs) > 1.0e-12)
    shs = -1.0 / shs;

  auto yyt = make_shared<Matrix>(size_, size_);
  auto sst = make_shared<Matrix>(size_, size_);
  dger_(size_, size_, shs, s->data(), 1, s->data(), 1, sst->data(), size_);
  dger_(size_, size_, ys, y->data(), 1, y->data(), 1, yyt->data(), size_);

  auto bsst = make_shared<Matrix>(*hess_ * *sst * *hess_);

  auto hess = make_shared<Matrix>(*bsst + *yyt);

  return hess;
}


shared_ptr<Matrix> Opt::hessian_update_psb(shared_ptr<const GradFile> y, shared_ptr<const GradFile> s, shared_ptr<const GradFile> z) const {
  // Hessian update with PSB
  double ss = s->dot_product(s);
  double ss2 = ss * ss;
  double sz = s->dot_product(z);

  ss = 1.0 / ss;
  ss2= -sz / ss2;

  auto szt = make_shared<Matrix>(size_, size_);
  auto zst = make_shared<Matrix>(size_, size_);
  auto sst = make_shared<Matrix>(size_, size_);
  dger_(size_, size_, ss, s->data(), 1, z->data(), 1, szt->data(), size_);
  dger_(size_, size_, ss, z->data(), 1, s->data(), 1, zst->data(), size_);
  dger_(size_, size_, ss2, s->data(), 1, s->data(), 1, sst->data(), size_);

  auto hess = make_shared<Matrix>(*szt + *zst + *sst);

  return hess;
}


shared_ptr<Matrix> Opt::hessian_update_bofill(shared_ptr<const GradFile> y, shared_ptr<const GradFile> s, shared_ptr<const GradFile> z) const {
  // Hessian update with bofill
  double ss = s->dot_product(s);
  double sz = s->dot_product(z);
  double zz = z->dot_product(z);

  double phi = sz * sz / (zz * ss);

  auto hess_psb = hessian_update_psb(y, s, z);
  auto hess_sr1 = hessian_update_sr1(y, s, z);
  hess_sr1->scale(phi);
  hess_psb->scale(1.0 - phi);

  auto hess = make_shared<Matrix>(*hess_psb + *hess_sr1);

  return hess;
}


tuple<double,double,shared_ptr<XYZFile>> Opt::get_step() const {
  auto displ = make_shared<XYZFile>(dispsize_);

  double predictedchange = 0.0, predictedchange_prev = 0.0;
  if (optinfo()->algorithm()->is_nr()) {
    displ = get_step_nr();
  } else if (optinfo()->algorithm()->is_rfo()) {
    tie(predictedchange, predictedchange_prev, displ) = get_step_rfo();
  } else if (optinfo()->algorithm()->is_ef()) {
    if (optinfo()->opttype()->is_transition())
      displ = get_step_ef_pn(0);
    else
      displ = get_step_ef();
  }
  return tie(predictedchange, predictedchange_prev, displ);
}


shared_ptr<XYZFile> Opt::get_step_nr() const {

  // Quasi-Newton-Raphson step
  auto displ = make_shared<XYZFile>(dispsize_);

  shared_ptr<Matrix> hinv(hess_);
  hinv->inverse();

  copy_n(grad_->data(), size_, displ->data());
  displ = displ->transform(hinv, /*transpose=*/false);

  displ->scale(-1.0);

  const double stepnorm = displ->norm();
  if (stepnorm > (maxstep_))
    displ->scale(maxstep_ / stepnorm);

  return displ;
}


shared_ptr<XYZFile> Opt::get_step_ef() const {
  // Eigenvector following by Baker
  auto displ = make_shared<XYZFile>(dispsize_);
  auto hess = make_shared<Matrix>(*hess_);

  // diagonalize Hessian
  VectorB eigv(size_);
  hess->diagonalize(eigv);

  // evaluate eigenvalue iteratively
  double lambda = 100.0;
  VectorB f(size_);

  for (int i = 0; i != size_; ++i) {
    copy_n(hess->element_ptr(0, i), size_, displ->data());
    f[i] = -displ->dot_product(grad_);
  }

  for (int iiter = 0; iiter != 100; ++iiter) {
    const double lambda_prev = lambda;
    double lambda_n = 0.0;
    for (int i = 0; i != size_; ++i)
      lambda_n += -(f[i] * f[i]) / (eigv[i] - lambda);
    lambda = lambda_n;

    const double error = fabs(lambda_prev - lambda);
    if (error < 1.0e-8)
      break;
  }

  displ->zero();
  for (int i = 0; i != size_; ++i) {
    auto dispb = make_shared<XYZFile>(dispsize_);
    const double fb = f[i] / (eigv[i] - lambda);
    copy_n(hess->element_ptr(0, i), size_, dispb->data());
    dispb->scale(fb);
    *displ += *dispb;
  }

  if (displ->norm() > maxstep_)
    displ->scale(maxstep_ / displ->norm());

  return displ;
}


shared_ptr<XYZFile> Opt::get_step_ef_pn(const int mode) const {
  // Eigenvector following by Baker
  auto displ = make_shared<XYZFile>(dispsize_);
  auto hess = make_shared<Matrix>(*hess_);

  // diagonalize Hessian
  VectorB eigv(size_);
  hess->diagonalize(eigv);

  cout << " Hessian Eigenvalues = " << endl;
  for (int i = 0; i != size_; ++i) {
    cout << setprecision(10) << eigv[i] << endl;
  }

  // partition lambda
  double lambda_p = 100.0;
  double lambda_n = 100.0;
  const int size_p = 1;
  const int size_n = size_ - size_p;

  VectorB f1(size_p);
  VectorB f2(size_n);

  copy_n(hess->element_ptr(0, mode), size_, displ->data());
  f1[0] = displ->dot_product(grad_);

  for (int j = 0, k = 0; j != size_; ++j) {
    if (j == mode) continue;
    copy_n(hess->element_ptr(0, j), size_, displ->data());
    f2[k] = displ->dot_product(grad_);
    ++k;
  }

  lambda_p = 0.5 * (eigv[mode] + sqrt(eigv[mode] * eigv[mode] + 4.0 * f1[0] * f1[0]));

  for (int iiter = 0; iiter != 500; ++iiter) {
    const double lambda_prev = lambda_n;
    double lambda_n_n = 0.0;
    for (int i = 0, k = 0; i != size_; ++i) {
      if (i == mode) continue;
      lambda_n_n += -(f2[k] * f2[k]) / (eigv[i] - lambda_n);
      ++k;
    }
    lambda_n = lambda_n_n;

    const double error = fabs(lambda_prev - lambda_n);
    cout << " error = " << setprecision(10) << error << setw(20) << lambda_n << endl;
    if (error < 1.0e-8)
      break;
  }

  displ->zero();

  for (int i = 0; i != size_p; ++i) {
    auto dispb = make_shared<XYZFile>(dispsize_);
    const double fb = -f1[i] / (eigv[mode] - lambda_p);
    copy_n(hess->element_ptr(0, mode), size_, dispb->data());
    dispb->scale(fb);
    *displ += *dispb;
  }

  for (int j = 0, k = 0; j != size_; ++j) {
    if (j == mode) continue;
    auto dispb = make_shared<XYZFile>(dispsize_);
    const double fb = -f2[k] / (eigv[j] - lambda_n);
    copy_n(hess->element_ptr(0, j), size_, dispb->data());
    dispb->scale(fb);
    *displ += *dispb;
    ++k;
  }

  if (displ->norm() > maxstep_)
    displ->scale(maxstep_ / displ->norm());

  return displ;
}


tuple<double,double,shared_ptr<XYZFile>> Opt::get_step_rfo() const {
  // Rational function optimization (aka augmented Hessian)
  // Here we do scale lambda to get step < steplength
  auto displ = make_shared<XYZFile>(dispsize_);
  {
    auto aughes = make_shared<Matrix>(size_ + 1,size_ + 1);
    VectorB eigv(size_+1);
    double lambda = 1.0;
    while (1) {
      aughes->zero();
      aughes->add_block(1.0 * lambda, 1, 1, size_, size_, hess_);
      aughes->add_block(1.0, 1, 0, size_, 1, grad_->data());
      aughes->add_block(1.0, 0, 1, 1, size_, grad_->data());

      aughes->diagonalize(eigv);
      aughes->scale(lambda / aughes->element(0, 0));

      if (!optinfo()->opttype()->is_transition())
        copy_n(aughes->element_ptr(1, 0), size_, displ->data());
      else
        copy_n(aughes->element_ptr(1, 1), size_, displ->data());

      if (displ->norm() < maxstep_)
        break;
      else
        lambda /= 1.2;
    }
  }

  double predictedchange = 0.0, predictedchange_prev = 0.0;
  if (optinfo()->adaptive()) {
    // When we use adaptive steplength, we should predict quadratic energy change

    const double qg  = displ->dot_product(grad_);
    auto hq = make_shared<GradFile>(*(displ->transform(hess_, /*transpose=*/false)));
    const double qhq = displ->dot_product(hq);
    predictedchange_prev = predictedchange_;
    predictedchange = qg + 0.5 * qhq;
  }

  return tie(predictedchange, predictedchange_prev, displ);
}


shared_ptr<XYZFile> Opt::iterate_displ() const {
  shared_ptr<Molecule> currentv = make_shared<Geometry>(*current_);
  auto displ = make_shared<XYZFile>(*displ_);
  auto dqc = make_shared<XYZFile>(*displ_);
  bool flag = false;
  shared_ptr<XYZFile> qc = displ->clone();
  if (optinfo()->redundant()) {
    copy_n(bmat_red_[2]->data(), size_, qc->data());
  } else {
    copy_n(current_->xyz()->data(), current_->natom()*3, qc->data());
    qc = qc->transform(bmat_[0], true);
  }
  cout << endl << "  === Displacement transformation iteration === " << endl << endl;
  Timer timer;

  array<shared_ptr<const Matrix>,3> bmat = bmat_;
  array<shared_ptr<const Matrix>,5> bmat_red = bmat_red_;
  for (int i = 0; i != optinfo()->maxiter(); ++i) {
    if (optinfo()->redundant()) {
      displ = displ->transform(bmat_red[1], false);
    } else {
      displ = displ->transform(bmat[1], false);
    }
    currentv = make_shared<Molecule>(*currentv, displ, true);
    if (optinfo()->redundant()) {
      vector<vector<int>> bondlist;
      tie(bondlist, bmat_red) = currentv->compute_redundant_coordinate(bondlist_);
    } else {
      bmat = currentv->compute_internal_coordinate(bmat[0], optinfo()->bonds(), optinfo()->opttype()->is_transition(), /*verbose=*/false);
    }
    auto qcurrent = make_shared<XYZFile>(dispsize_);
    if (optinfo()->redundant()) {
      copy_n(bmat_red[2]->data(), size_, qcurrent->data());
      for (int i = 0; i != size_; ++i) {
        const double thresh = pi__ - 0.1;
        if (qcurrent->element(i%3, i/3) > thresh && qc->element(i%3, i/3) < -thresh)      qcurrent->element(i%3, i/3) -= 2.0 * pi__;
        else if (qcurrent->element(i%3, i/3) < -thresh && qc->element(i%3, i/3) > thresh) qcurrent->element(i%3, i/3) += 2.0 * pi__;
      }
    } else {
      copy_n(currentv->xyz()->data(), currentv->natom()*3, qcurrent->data());
      qcurrent = qcurrent->transform(bmat[0], true);
    }
    auto qdiff = make_shared<XYZFile>(*qcurrent - *qc);
    displ = make_shared<XYZFile>(*dqc - *qdiff);
    if ((displ->norm() / size_ > 0.1 || i == (optinfo()->maxiter() - 1))) {
      cout << "  * Calculated displacement too large, switching to first-order" << endl;
      flag = true;
      break;
    }
    cout << setw(7) << i << setprecision(10) << setw(15) << displ->norm() / (currentv->natom()*3) << setw(10) << setprecision(2) << timer.tick() << endl;
    if (displ->norm() < 1.0e-7)
      break;
  }

  if (flag) {
    if (optinfo()->redundant()) {
      displ = displ_->transform(bmat_red[1], false);
    } else {
      displ = displ_->transform(bmat[1], false);
    }
  } else {
    *displ = *(currentv->xyz()) - *(current_->xyz());
  }
  cout << endl << endl;
  return displ;
}


double Opt::do_adaptive(const int iter) const {
  // Fletcher's adaptive stepsize algorithm, works with rfo

  const bool flag = iter > 1 && optinfo()->algorithm()->is_rfo() && optinfo()->opttype()->is_energy();

  double maxstep = maxstep_;
  if (flag) {
    const double realchange = en_ - prev_en_.back();
    double predreal_ratio = realchange / predictedchange_prev_;
    if (predreal_ratio > 1.0)
      predreal_ratio = 1.0 / predreal_ratio;

    if (predreal_ratio > 0.75 && (displ_->norm() > (0.80 * maxstep)))
      maxstep *= 2.0;
    else if (predreal_ratio < 0.25)
      maxstep *= 0.25;
  }

  return maxstep;
}

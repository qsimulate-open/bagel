//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: get_grad.cc
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
#include <src/wfn/get_energy.h>
#include <src/opt/opt.h>
#include <src/grad/finite.h>

using namespace std;
using namespace bagel;

tuple<double,double,shared_ptr<const Reference>,shared_ptr<GradFile>> Opt::get_mecigrad(shared_ptr<PTree> cinput, shared_ptr<const Reference> ref) const {
  // MECI: minimize E2 while E2-E1 = 0 -> project out g and h from D(E2). [Bearpark, Robb, Schlegel (CPL 1994, 223, 269)]

  auto out = make_shared<GradFile>(current_->natom());
  int n3 = current_->natom() * 3;
  double en1 = 0.0, en2 = 0.0;
  shared_ptr<GradFile> cgrad1;
  shared_ptr<GradFile> cgrad2;
  shared_ptr<GradFile> x2;

  shared_ptr<const Reference> prev_ref;
  if (method_ == "casscf") {
    GradEval<CASSCF> eval1(cinput, current_, ref);
    cgrad1 = make_shared<GradFile>(*eval1.compute("force", make_shared<GradInfo>(optinfo(), optinfo()->target_state())));
    prev_ref = eval1.ref();
    en2 = eval1.energy();

    cgrad2 = make_shared<GradFile>(*eval1.compute("force", make_shared<GradInfo>(optinfo(), optinfo()->target_state2())));
    en1 = eval1.energy();

    x2 = make_shared<GradFile>(*eval1.compute("nacme", optinfo()));
  } else if (method_ == "caspt2") {
    GradEval<CASPT2Grad> eval1(cinput, current_, ref);
    cgrad1 = make_shared<GradFile>(*eval1.compute("force", make_shared<GradInfo>(optinfo(), optinfo()->target_state())));
    prev_ref = eval1.ref();
    en2 = eval1.energy();

    cgrad2 = make_shared<GradFile>(*eval1.compute("force", make_shared<GradInfo>(optinfo(), optinfo()->target_state2())));
    en1 = eval1.energy();

    x2 = make_shared<GradFile>(*eval1.compute("nacme", optinfo()));
  } else {
    throw logic_error ("Conical intersection search currently only available for CASSCF or CASPT2");
  }

  if (optinfo()->qmmm()) {
    double mmen;
    shared_ptr<GradFile> mmgrad;
    qmmm_driver_->edit_input(current_);
    tie(mmen, mmgrad) = qmmm_driver_->do_grad(current_->natom());
    *cgrad2 = *cgrad2 + *mmgrad;
    *cgrad1 = *cgrad1 + *mmgrad;
    en1 += mmen;
    en2 += mmen;
  }

  auto x1 = make_shared<GradFile>(*cgrad1 - *cgrad2);
  const double x1norm = x1->norm();
  if (x1norm > 1.0e-8)
    x1->scale(1.0 / x1norm);
  auto xf = make_shared<GradFile>(*x1);
  const double en  = en2 - en1;
  xf->scale(2.0 * en / x1norm);
  auto xg = make_shared<GradFile>(*cgrad1);
  {
    const double x2norm = x2->norm();
    if (x2norm > 1.0e-8)
      x2->scale(1.0 / x2norm);
  }

  auto proj = make_shared<Matrix>(n3, n3);
  proj->unit();
  dger_(n3, n3, -1.0, x1->data(), 1, x1->data(), 1, proj->data(), n3);
  x2 = x2->transform(proj, false);
  {
    const double x2norm = x2->norm();
    if (x2norm > 1.0e-8)
      x2->scale(1.0 / x2norm);
  }
  proj->unit();
  dger_(n3, n3, -1.0, x1->data(), 1, x1->data(), 1, proj->data(), n3);
  dger_(n3, n3, -1.0, x2->data(), 1, x2->data(), 1, proj->data(), n3);
  xg = xg->transform(proj, /*transpose=*/false);
  *out = optinfo()->thielc3() * (*xf * optinfo()->thielc4() + *xg * (1.0 - optinfo()->thielc4()));

  return tie(en2, en, prev_ref, out);
}


tuple<double,shared_ptr<GradFile>> Opt::get_euclidean_dist(shared_ptr<const XYZFile> a, shared_ptr<const XYZFile> ref) const {
  // This aligns two structures and evaluates q^2 and dq^2 / dX [Rhee (JCP 2000, 113, 6021)].
  const int natom = current_->natom();
  auto q_eckt = make_shared<XYZFile>(*ref);

  // translate cog of q_eckt into a
  vector<double> ref_center;
  vector<double> a_center;

  ref_center.assign(3, 0.0);
  a_center.assign(3, 0.0);
  for (int iatom = 0; iatom != natom; ++iatom) {
    for (int j = 0; j != 3; ++j) {
      a_center[j] += a->element(j,iatom) / natom;
      ref_center[j] += q_eckt->element(j,iatom) / natom;
    }
  }
  ref_center[0] -= a_center[0];
  ref_center[1] -= a_center[1];
  ref_center[2] -= a_center[2];

  for (int iatom = 0; iatom != natom; ++iatom)
    for (int j = 0; j != 3; ++j)
      q_eckt->element(j,iatom) -= ref_center[j];

  // make q_eckt rotated to match with current geometry
  for (int iiter = 0; iiter != 500; ++iiter) {
    auto fmat = make_shared<Matrix>(3,3);
    fmat->fill(0.0);
    for (int iatom = 0; iatom != natom; ++iatom) {
      fmat->element(0,0) += q_eckt->element(0,iatom) * a->element(0,iatom);
      fmat->element(1,0) += q_eckt->element(0,iatom) * a->element(1,iatom);
      fmat->element(2,0) += q_eckt->element(0,iatom) * a->element(2,iatom);
      fmat->element(0,1) += q_eckt->element(1,iatom) * a->element(0,iatom);
      fmat->element(1,1) += q_eckt->element(1,iatom) * a->element(1,iatom);
      fmat->element(2,1) += q_eckt->element(1,iatom) * a->element(2,iatom);
      fmat->element(0,2) += q_eckt->element(2,iatom) * a->element(0,iatom);
      fmat->element(1,2) += q_eckt->element(2,iatom) * a->element(1,iatom);
      fmat->element(2,2) += q_eckt->element(2,iatom) * a->element(2,iatom);
    }
    const vector<double> grad { fmat->element(2,1) - fmat->element(1,2),
                                fmat->element(0,2) - fmat->element(2,0),
                                fmat->element(1,0) - fmat->element(0,1) };
    const double hdiag0 = -(fmat->element(1,1) + fmat->element(2,2));
    const double hdiag1 = -(fmat->element(2,2) + fmat->element(0,0));
    const double hdiag2 = -(fmat->element(0,0) + fmat->element(1,1));

    const double deth = hdiag0 * hdiag1 * hdiag2 + 2.0 * fmat->element(1,0)*fmat->element(2,1)*fmat->element(2,0)
                       -hdiag0 * fmat->element(2,1) * fmat->element(2,1) - hdiag1 * fmat->element(2,0) * fmat->element(2,0)
                       -hdiag2 * fmat->element(1,0) * fmat->element(1,0);

    const double norm = grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2];
    auto tfm = make_shared<Matrix>(3,3);
    if (norm < 1.0e-8) {
      const double trac = hdiag0 + hdiag1 + hdiag2;
      const double det3 = (hdiag0+hdiag2)*hdiag1 + hdiag0*hdiag2 - fmat->element(1,0)*fmat->element(1,0) - fmat->element(2,1)*fmat->element(2,1) - fmat->element(2,0)*fmat->element(2,0);
      if (deth < 1.0e-8 && trac < 0.0 && det3 > 0.0)
        break;         // reached convergence
      // reached local minima
      fmat->element(0,0) = hdiag0;
      fmat->element(1,1) = hdiag1;
      fmat->element(2,2) = hdiag2;
      auto eig = make_shared<VectorB>(3);
      fmat->diagonalize(*eig);
      for (int i = 0; i != 3; ++i) {
        for (int j = 0; j != 3; ++j)
          tfm->element(i,j) = 2.0 * fmat->element(i,2) * fmat->element(j,2);
        tfm->element(i,i) -= 1.0;
      }
    } else {
      auto hinv = make_shared<Matrix>(3,3);
      double alpha, beta, gamma;
      if (fabs(deth) > 1.0e-8) {
        hinv->element(0,0) = hdiag1 * hdiag2 - fmat->element(2,1) * fmat->element(2,1);
        hinv->element(1,1) = hdiag0 * hdiag2 - fmat->element(2,0) * fmat->element(2,0);
        hinv->element(2,2) = hdiag0 * hdiag1 - fmat->element(1,0) * fmat->element(1,0);
        hinv->element(1,0) = fmat->element(2,0) * fmat->element(2,1) - fmat->element(1,0) * hdiag2;
        hinv->element(2,0) = fmat->element(1,0) * fmat->element(2,1) - fmat->element(2,0) * hdiag1;
        hinv->element(2,1) = fmat->element(2,0) * fmat->element(1,0) - fmat->element(2,1) * hdiag0;
        alpha = -(hinv->element(0,0) * grad[0] + hinv->element(1,0) * grad[1] + hinv->element(2,0) * grad[2]) / deth;
        beta  = -(hinv->element(1,0) * grad[0] + hinv->element(1,1) * grad[1] + hinv->element(2,1) * grad[2]) / deth;
        gamma = -(hinv->element(2,0) * grad[0] + hinv->element(2,1) * grad[1] + hinv->element(2,2) * grad[2]) / deth;
      } else {
        fmat->element(0,0) = hdiag0;
        fmat->element(1,1) = hdiag1;
        fmat->element(2,2) = hdiag2;
        fmat->element(1,2) = fmat->element(2,1);
        fmat->element(0,2) = fmat->element(2,0);
        fmat->element(0,1) = fmat->element(1,0);
        auto eig = make_shared<VectorB>(3);
        fmat->diagonalize(*eig);
        for (int i = 0; i != 3; ++i)
          if (fabs((*eig)[i]) < 1.0e-8)
            (*eig)[i] = 0.0;
          else
            (*eig)[i] = 1.0 / (*eig)[i];
        for (int i = 0; i != 3; ++i)
          for (int j = 0; j <= i; j++)
            hinv->element(i,j) = fmat->element(i,0) * (*eig)[0] * fmat->element(j,0)
                               + fmat->element(i,1) * (*eig)[1] * fmat->element(j,1)
                               + fmat->element(i,2) * (*eig)[2] * fmat->element(j,2);
        alpha = -(hinv->element(0,0) * grad[0] + hinv->element(1,0) * grad[1] + hinv->element(2,0) * grad[2]);
        beta  = -(hinv->element(1,0) * grad[0] + hinv->element(1,1) * grad[1] + hinv->element(2,1) * grad[2]);
        gamma = -(hinv->element(2,0) * grad[0] + hinv->element(2,1) * grad[1] + hinv->element(2,2) * grad[2]);
      }

      const double cosa = cos(alpha);
      const double sina = sin(alpha);
      const double cosb = cos(beta);
      const double sinb = sin(beta);
      const double cosc = cos(gamma);
      const double sinc = sin(gamma);
      tfm->element(0,0) =  cosb * cosc;
      tfm->element(1,0) =  sina * sinb * cosc + cosa * sinc;
      tfm->element(2,0) =  sina * sinc - cosa * sinb * cosc;
      tfm->element(0,1) = -cosb * sinc;
      tfm->element(1,1) =  cosa * cosc - sina * sinb * sinc;
      tfm->element(2,1) =  cosa * sinb * sinc + sina * cosc;
      tfm->element(0,2) =  sinb;
      tfm->element(1,2) = -sina * cosb;
      tfm->element(2,2) =  cosa * cosb;
    }
    for (int iatom = 0; iatom != natom; ++iatom) {
      vector<double> qdum(3);
      for (int j = 0; j != 3; ++j)
        qdum[j] = tfm->element(j,0) * q_eckt->element(0,iatom) + tfm->element(j,1) * q_eckt->element(1,iatom) + tfm->element(j,2) * q_eckt->element(2,iatom);
      for (int j = 0; j != 3; ++j)
        q_eckt->element(j,iatom) = qdum[j];
    }
  }

  // q^2 = sum_{iatom} [a - q_eckt]^2
  double q2 = 0.0;
  for (int iatom = 0; iatom != natom; ++iatom)
    for (int j = 0; j != 3; ++j)
      q2 += (a->element(j,iatom) - q_eckt->element(j,iatom)) * (a->element(j,iatom) - q_eckt->element(j,iatom));

  // dq^2/dx = 2 * [a - q_eckt]
  auto dqdx = make_shared<GradFile>(natom);
  for (int iatom = 0; iatom != natom; ++iatom)
    for (int j = 0; j != 3; ++j)
      dqdx->element(j,iatom) = 2.0 * (a->element(j,iatom) - q_eckt->element(j,iatom));

  double dist;
  if (q2 < 1.0e-12)
    dist = 0.0;
  else
    dist = sqrt(q2);

  // Scale dq/dx by (1.0 / (dist + 0.1) - 0.5 * dist / (dist + 0.1)^2) [differentiation of the first term in eq. 15 in JPCB 2008, 112, 405]
  if (dist > 1.0e-6)
    dqdx->scale(1.0 / (dist + 0.1/au2angstrom__) - 0.5 * dist / ((dist + 0.1/au2angstrom__) * (dist + 0.1/au2angstrom__)));

  return tie(dist, dqdx);
}


tuple<double,double,shared_ptr<const Reference>,shared_ptr<GradFile>> Opt::get_mdcigrad(shared_ptr<PTree> cinput, shared_ptr<const Reference> ref) const {
  // Concept of MDCI is suggested by Levine, Coe and Martinez (JPCB 2008, 112, 405) -- CI geometry with minimized distance to reference geometry.
  // Original article used "target function" to minimize.
  // Here, we make use of gradient projection :: minimize distance while E2-E1 = 0 -> project out g and h from the vector d(distance).

  auto out = make_shared<GradFile>(current_->natom());
  int n3 = current_->natom() * 3;
  double en1 = 0.0, en2 = 0.0;
  shared_ptr<GradFile> cgrad1;
  shared_ptr<GradFile> cgrad2;
  shared_ptr<GradFile> x2;

  shared_ptr<const Reference> prev_ref;
  if (method_ == "casscf") {
    GradEval<CASSCF> eval1(cinput, current_, ref);
    cgrad1 = make_shared<GradFile>(*eval1.compute("force", make_shared<GradInfo>(optinfo(), optinfo()->target_state())));
    prev_ref = eval1.ref();
    en2 = eval1.energy();

    cgrad2 = make_shared<GradFile>(*eval1.compute("force", make_shared<GradInfo>(optinfo(), optinfo()->target_state2())));
    en1 = eval1.energy();

    x2 = make_shared<GradFile>(*eval1.compute("nacme", optinfo()));
  } else if (method_ == "caspt2") {
    GradEval<CASPT2Grad> eval1(cinput, current_, ref);
    cgrad1 = make_shared<GradFile>(*eval1.compute("force", make_shared<GradInfo>(optinfo(), optinfo()->target_state())));
    prev_ref = eval1.ref();
    en2 = eval1.energy();

    cgrad2 = make_shared<GradFile>(*eval1.compute("force", make_shared<GradInfo>(optinfo(), optinfo()->target_state2())));
    en1 = eval1.energy();

    x2 = make_shared<GradFile>(*eval1.compute("nacme", optinfo()));
  } else {
    throw logic_error ("Conical intersection search currently only available for CASSCF or CASPT2");
  }

  if (optinfo()->qmmm()) {
    double mmen;
    shared_ptr<GradFile> mmgrad;
    qmmm_driver_->edit_input(current_);
    tie(mmen, mmgrad) = qmmm_driver_->do_grad(current_->natom());
    *cgrad2 = *cgrad2 + *mmgrad;
    *cgrad1 = *cgrad1 + *mmgrad;
    en1 += mmen;
    en2 += mmen;
  }

  auto x1 = make_shared<GradFile>(*cgrad1 - *cgrad2);
  const double x1norm = x1->norm();
  if (x1norm > 1.0e-8)
    x1->scale(1.0 / x1norm);
  auto xf = make_shared<GradFile>(*x1);
  const double en  = en2 - en1;
  xf->scale(2.0 * en / x1norm);

  shared_ptr<GradFile> xg;
  const bool refg = idata_->get<bool>("mdci_reference_geometry", false);
  shared_ptr<XYZFile> ref_xyz;
  if (refg) {
    auto input = idata_->get_child("refgeom");
    auto geom_ref = make_shared<Geometry>(input);
    ref_xyz = make_shared<XYZFile>(*(geom_ref->xyz()));
  } else {
    ref_xyz = make_shared<XYZFile>(*prev_xyz_[0]);
  }

  double dist;
  tie(dist, xg) = get_euclidean_dist(current_->xyz(), ref_xyz);
  {
    const double x2norm = x2->norm();
    if (x2norm > 1.0e-8)
      x2->scale(1.0 / x2norm);
  }

  auto proj = make_shared<Matrix>(n3, n3);
  proj->unit();
  dger_(n3, n3, -1.0, x1->data(), 1, x1->data(), 1, proj->data(), n3);
  x2 = x2->transform(proj, false);
  {
    const double x2norm = x2->norm();
    if (x2norm > 1.0e-8)
      x2->scale(1.0 / x2norm);
  }
  proj->unit();
  dger_(n3, n3, -1.0, x1->data(), 1, x1->data(), 1, proj->data(), n3);
  dger_(n3, n3, -1.0, x2->data(), 1, x2->data(), 1, proj->data(), n3);
  xg = xg->transform(proj, /*transpose=*/false);
  *out = optinfo()->thielc3() * (*xf * optinfo()->thielc4() + *xg * (1.0 - optinfo()->thielc4()));

  return tie(dist, en, prev_ref, out);
}


tuple<double,shared_ptr<const Reference>,shared_ptr<GradFile>> Opt::get_grad_energy(shared_ptr<PTree> cinput, shared_ptr<const Reference> ref) const {
  auto out = make_shared<GradFile>(current_->natom());
  shared_ptr<const Reference> prev_ref;
  double en;
  bool numerical = optinfo()->numerical();

  if (!numerical) {
    if (method_ == "uhf") {

      GradEval<UHF> eval(cinput, current_, ref);
      out = eval.compute("force", optinfo());
      prev_ref = eval.ref();
      en = eval.energy();

    } else if (method_ == "rohf") {

      GradEval<ROHF> eval(cinput, current_, ref);
      out = eval.compute("force", optinfo());
      prev_ref = eval.ref();
      en = eval.energy();

    } else if (method_ == "hf") {

      GradEval<RHF> eval(cinput, current_, ref);
      out = eval.compute("force", optinfo());
      prev_ref = eval.ref();
      en = eval.energy();

    } else if (method_ == "ks") {

      GradEval<KS> eval(cinput, current_, ref);
      out = eval.compute("force", optinfo());
      prev_ref = eval.ref();
      en = eval.energy();

    } else if (method_ == "dhf") {

      GradEval<Dirac> eval(cinput, current_, ref);
      out = eval.compute("force", optinfo());
      prev_ref = eval.ref();
      en = eval.energy();

    } else if (method_ == "mp2") {

      GradEval<MP2Grad> eval(cinput, current_, ref);
      out = eval.compute("force", optinfo());
      prev_ref = eval.ref();
      en = eval.energy();

    } else if (method_ == "casscf") {

      GradEval<CASSCF> eval(cinput, current_, ref);
      out = eval.compute("force", optinfo());
      prev_ref = eval.ref();
      en = eval.energy();

    } else if (method_ == "caspt2") {

      GradEval<CASPT2Grad> eval(cinput, current_, ref);
      out = eval.compute("force", optinfo());
      prev_ref = eval.ref();
      en = eval.energy();

    } else {

      cout << "   There is no analytical gradient available. Numerical gradient will be used." << endl;
      numerical = true;

    }
  }

  if (numerical) {
    auto m = idata_->get_child("method");
    const int nproc = idata_->get<int>("nproc", 1);
    const double dx = idata_->get<double>("numerical_dx", 0.001);
    FiniteGrad eval(m, current_, ref, optinfo()->target_state(), dx, nproc);
    out = eval.compute();
    prev_ref = eval.ref();
    en = eval.energy();
  }

  if (optinfo()->qmmm()) {
    double mmen;
    shared_ptr<GradFile> mmgrad;
    qmmm_driver_->edit_input(current_);
    tie(mmen, mmgrad) = qmmm_driver_->do_grad(current_->natom());
    *out = *out + *mmgrad;
    en += mmen;
  }

  return tie(en, prev_ref, out);
}


tuple<double,double,shared_ptr<const Reference>,shared_ptr<GradFile>> Opt::get_grad(shared_ptr<PTree> cinput, shared_ptr<const Reference> ref) const {
  auto out = make_shared<GradFile>(current_->natom());

  shared_ptr<const Reference> prev_ref;
  double param1,param2;
  if (optinfo()->opttype()->is_mdci()) {
    tie(param1, param2, prev_ref, out) = get_mdcigrad(cinput, ref);
  } else if (optinfo()->opttype()->is_conical()) {
    tie(param1, param2, prev_ref, out) = get_mecigrad(cinput, ref);
  } else {
    tie(param1, prev_ref, out) = get_grad_energy(cinput, ref);
  }

  return tie(param1, param2, prev_ref, out);
}


tuple<shared_ptr<PTree>,shared_ptr<const Reference>,shared_ptr<const Geometry>> Opt::get_grad_input() const {
  shared_ptr<PTree> cinput;
  shared_ptr<const Reference> ref;
  auto current = make_shared<const Geometry>(*current_);

  if (prev_ref_ && prev_ref_->is_skelton())
    throw runtime_error("Perform a single-point calculation before optimization when starting from Molden");

  if (!prev_ref_ || optinfo()->scratch()) {
    auto m = input_->begin();
    for ( ; m != --input_->end(); ++m) {
      const string title = to_lower((*m)->get<string>("title", ""));
      if (title != "molecule") {
        tie(ignore, ref) = get_energy(title, *m, current, ref);
      } else {
        current = make_shared<const Geometry>(*current, *m);
        if (ref)
          ref = ref->project_coeff(current);
      }
    }
    cinput = make_shared<PTree>(**m);
  } else {
    ref = prev_ref_->project_coeff(current);
    cinput = make_shared<PTree>(**input_->rbegin());
  }
  cinput->put("_gradient", true);

  return tie(cinput, ref, current);
}

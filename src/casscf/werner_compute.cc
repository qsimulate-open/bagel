//
// BAGEL - Parallel electron correlation program.
// Filename: werner_compute.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/casscf/werner.h>
#include <src/casscf/jkop.h>

using namespace std;
using namespace bagel;

#define DEBUG4INDEX



shared_ptr<Matrix> WernerKnowles::compute_sigma_R(const shared_ptr<const Jvec> jvec, const shared_ptr<const Matrix> dR,
                const shared_ptr<const Matrix> C, const shared_ptr<const Matrix> U) {

  // compute U dR
  shared_ptr<Matrix> UdR(new Matrix(*U**dR));

  // compute  1/2 (C dR + dR C)
  shared_ptr<Matrix> dRA(new Matrix(*C**dR+*dR**C));

  // update B
//#ifndef DEBUG4INDEX
#if 1
  shared_ptr<Matrix> new_bvec = compute_bvec(jvec, UdR, UdR, coeff_);
#else
  shared_ptr<Matrix> gg(new Matrix(*coeff_ * *UdR));
  shared_ptr<Matrix> tmp = jk.contract(gg);
  shared_ptr<Matrix> new_bvec(new Matrix(*coeff_ % *tmp));
#endif

  // compute U^dagger B - B^dagger U
  // compute Eq.29
  shared_ptr<Matrix> sigma(new Matrix(*U%*new_bvec-*new_bvec%*U-*dRA));
  sigma->purify_redrotation(nclosed_,nact_,nvirt_);

  return sigma;
}



shared_ptr<const Matrix> WernerKnowles::compute_denom(const shared_ptr<const Matrix> C) {
#ifndef DEBUG4INDEX
  shared_ptr<Matrix> f;
  shared_ptr<QFile>    fact, factp, gaa;
  shared_ptr<RotFile> denom_;
  one_body_operators(f, fact, factp, gaa, denom_, false);
  shared_ptr<Matrix> denom = denom_->unpack_sym(geom_, 1.0e1);
#else
  JKop jk(geom_->df(), coeff_, hcore_, fci_, nocc_, nclosed_, nact_);
  unique_ptr<double[]> cdiag = C->diag();
  shared_ptr<Matrix> denom = jk.denom();
  for (int i = 0; i != nocc_; ++i) {
    for (int j = nocc_; j != nbasis_; ++j) {
      denom->element(j, i) -= cdiag[i] + cdiag[j];
      denom->element(i, j) -= cdiag[i] + cdiag[j];
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = nclosed_; j != nocc_; ++j) {
      denom->element(j, i) -= cdiag[i] + cdiag[j];
      denom->element(i, j) -= cdiag[i] + cdiag[j];
    }
  }
#endif
  return denom;
}


// compute B according to Eq. (19).
// B = 2[h_rs U_sj D_ji + (rs|D)Usj <D|ji> + 2(rk|D)(D|ls)T_sj D_jl,ik]
shared_ptr<Matrix> WernerKnowles::compute_bvec(const shared_ptr<const Jvec> jvec, shared_ptr<Matrix> u, const shared_ptr<const Coeff> cc) {
  shared_ptr<Matrix> t(new Matrix(*u));
  for (int i = 0; i != u->ndim(); ++i) t->element(i,i) -= 1.0;
  return compute_bvec(jvec, u, t, cc);
}




shared_ptr<Matrix> WernerKnowles::compute_bvec(const shared_ptr<const Jvec> jvec,
                                                 shared_ptr<Matrix> u, shared_ptr<Matrix> t, const shared_ptr<const Coeff> cc) {
  shared_ptr<const DFDist> df = geom_->df();
  const int naux = df->naux();
  const int nbas = df->nbasis0();
  assert(df->nbasis0() == df->nbasis1());

// TODO make sure if this works
if (nbasis_ != nbas) throw runtime_error("I should examine this case...");

  shared_ptr<Matrix> out(new Matrix(geom_->nbasis(), geom_->nbasis()));
  {
    shared_ptr<Matrix> hcore_mo_(new Matrix(*cc % *hcore_ * *cc));

    // first term
    Matrix all1(geom_->nbasis(), geom_->nbasis());
    for (int i = 0; i != nclosed_; ++i) all1.element(i,i) = 2.0;
    for (int i = 0; i != nact_; ++i) copy_n(fci_->rdm1_av()->data()+nact_*i, nact_, all1.element_ptr(nclosed_, i+nclosed_));
    shared_ptr<Matrix> buf(new Matrix(geom_->nbasis(), geom_->nbasis()));
    dgemm_("N", "N", nbasis_, nocc_, nbasis_, 1.0, hcore_mo_->data(), nbasis_, u->data(), nbas, 0.0, buf->data(), nbas);
    dgemm_("N", "N", nbasis_, nocc_, nocc_, 2.0, buf->data(), nbas, all1.data(), nbas, 0.0, out->data(), nbas);
  }

  unique_ptr<double[]> tmp;
  // second term
  {
    Matrix Umat(*cc * *u);
    shared_ptr<DFHalfDist> half = df->compute_half_transform(Umat.data(), nocc_);
    tmp = move(half->form_2index(jvec->jvec(), 2.0));
  }

  // third term
  if (t->norm() > 1.0e-15) {
    Matrix Tmat(*cc * *t);
    shared_ptr<DFFullDist> full = jvec->half()->compute_second_transform(Tmat.data(), nocc_)->apply_2rdm(jvec->rdm2_all());
    unique_ptr<double[]> tmp2 = jvec->half()->form_2index(full, 4.0);
    daxpy_(nbas*nocc_, 1.0, tmp2, 1, tmp, 1);
  }

  dgemm_("T", "N", nbasis_, nocc_, nbas, 1.0, cc->data(), nbas, tmp.get(), nbas, 1.0, out->data(), nbasis_);

  return out;
}



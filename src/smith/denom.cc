//
// BAGEL - Parallel electron correlation program.
// Filename: denom.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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


#include <src/smith/denom.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;
using namespace SMITH;

Denom::Denom(shared_ptr<const Matrix> fock, const int nstates, const double th = 1.0e-8) : fock_(fock), thresh_(th) {
  const size_t ndim = fock->mdim() * nstates;
  const size_t ndim2 = ndim * ndim;
  const size_t ndim3 = ndim2 * ndim;
  assert(ndim == fock->ndim());

  shalf_x_ = make_shared<Matrix>(ndim, ndim);
  shalf_h_ = make_shared<Matrix>(ndim, ndim);
  shalf_xx_ = make_shared<Matrix>(ndim2, ndim2);
  shalf_hh_ = make_shared<Matrix>(ndim2, ndim2);
  shalf_xh_ = make_shared<Matrix>(ndim2*2, ndim2*2);
  shalf_xhh_ = make_shared<Matrix>(ndim3, ndim3);
  shalf_xxh_ = make_shared<Matrix>(ndim3, ndim3);

  work_x_ = make_shared<Matrix>(ndim, ndim);
  work_h_ = make_shared<Matrix>(ndim, ndim);
  work_xx_ = make_shared<Matrix>(ndim2, ndim2);
  work_hh_ = make_shared<Matrix>(ndim2, ndim2);
  work_xh_ = make_shared<Matrix>(ndim2*2, ndim2*2);
  work_xhh_ = make_shared<Matrix>(ndim3, ndim3);
  work_xxh_ = make_shared<Matrix>(ndim3, ndim3);

  denom_x_ = VectorB(ndim);
  denom_h_ = VectorB(ndim);
  denom_xx_ = VectorB(ndim2);
  denom_hh_ = VectorB(ndim2);
  denom_xh_ = VectorB(ndim2*2);
  denom_xhh_ = VectorB(ndim3);
  denom_xxh_ = VectorB(ndim3);
}


void Denom::append(const int jst, const int ist, shared_ptr<const RDM<1>> rdm1, shared_ptr<const RDM<2>> rdm2, shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4) {
  init_x_(jst, ist, rdm1, rdm2, rdm3, rdm4);
  init_h_(jst, ist, rdm1, rdm2, rdm3, rdm4);
  init_xx_(jst, ist, rdm1, rdm2, rdm3, rdm4);
  init_hh_(jst, ist, rdm1, rdm2, rdm3, rdm4);
  init_xh_(jst, ist, rdm1, rdm2, rdm3, rdm4);
  init_xhh_(jst, ist, rdm1, rdm2, rdm3, rdm4);
  init_xxh_(jst, ist, rdm1, rdm2, rdm3, rdm4);
}


void Denom::compute() {
  {
    shalf_x_->inverse_half(thresh_);
    Matrix tmp(*shalf_x_ % *work_x_ * *shalf_x_);
    tmp.diagonalize(denom_x_);
    shalf_x_ = make_shared<Matrix>(tmp % *shalf_x_);
  }
  {
    shalf_h_->inverse_half(thresh_);
    Matrix tmp(*shalf_h_ % *work_h_ * *shalf_h_);
    tmp.diagonalize(denom_h_);
    shalf_h_ = make_shared<Matrix>(tmp % *shalf_h_);
  }
  {
    shalf_xx_->inverse_half(thresh_);
    Matrix tmp(*shalf_xx_ % *work_xx_ * *shalf_xx_);
    tmp.diagonalize(denom_xx_);
    shalf_xx_ = make_shared<Matrix>(tmp % *shalf_xx_);
  }
  {
    shalf_hh_->inverse_half(thresh_);
    Matrix tmp(*shalf_hh_ % *work_hh_ * *shalf_hh_);
    tmp.diagonalize(denom_hh_);
    shalf_hh_ = make_shared<Matrix>(tmp % *shalf_hh_);
  }
  {
    shalf_xh_->inverse_half(thresh_);
    Matrix tmp(*shalf_xh_ % *work_xh_ * *shalf_xh_);
    tmp.diagonalize(denom_xh_);
    shalf_xh_ = make_shared<Matrix>(tmp % *shalf_xh_);
  }
  {
    shalf_xhh_->inverse_half(thresh_);
    Matrix tmp(*shalf_xhh_ % *work_xhh_ * *shalf_xhh_);
    tmp.diagonalize(denom_xhh_);
    shalf_xhh_ = make_shared<Matrix>(tmp % *shalf_xhh_);
  }
  {
    shalf_xxh_->inverse_half(thresh_);
    Matrix tmp(*shalf_xxh_ % *work_xxh_ * *shalf_xxh_);
    tmp.diagonalize(denom_xxh_);
    shalf_xxh_ = make_shared<Matrix>(tmp % *shalf_xxh_);
  }
}


void Denom::init_x_(const int jst, const int ist, shared_ptr<const RDM<1>> rdm1, shared_ptr<const RDM<2>> rdm2, shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4) {
  const size_t nact = rdm1->norb();
  const size_t dim = nact;
  const size_t size = dim*dim;

  shalf_x_->copy_block(dim*jst, dim*ist, dim, dim, rdm1);

  Matrix work2(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, rdm2->data(), size, fock_->data(), 1, 0.0, work2.data(), 1);
  work_x_->copy_block(dim*jst, dim*ist, dim, dim, work2);
}


void Denom::init_h_(const int jst, const int ist, shared_ptr<const RDM<1>> rdm1, shared_ptr<const RDM<2>> rdm2, shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact;
  const size_t size = dim*dim;
  auto shalf = make_shared<Matrix>(dim, dim);
  copy_n(rdm1->data(), size, shalf->data());
  shalf->scale(-1.0);
  shalf->add_diag(2.0); //.. making hole 1RDM
  shalf_h_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  shared_ptr<RDM<2>> ovl = rdm2->copy();
  ovl->scale(-1.0);
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k) {
        // see Celani eq. A7
        ovl->element(k, i, i, j) +=       shalf->element(k,j);
        ovl->element(i, k, j, i) += -1.0 * rdm1->element(k,j);
        ovl->element(i, i, k, j) +=  2.0 * rdm1->element(k,j);
      }
  Matrix work2(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, ovl->data(), size, fock_->data(), 1, 0.0, work2.data(), 1);
  work_h_->copy_block(dim*jst, dim*ist, dim, dim, work2);
}



void Denom::init_xx_(const int jst, const int ist, shared_ptr<const RDM<1>> rdm1, shared_ptr<const RDM<2>> rdm2, shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact;
  const size_t size = dim*dim;
  auto shalf = make_shared<Matrix>(dim, dim);
  shared_ptr<RDM<2>> tmp = rdm2->copy();
  sort_indices<0,2,1,3,0,1,1,1>(tmp->data(), shalf->data(), nact, nact, nact, nact);
  shalf_xx_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  auto work2 = make_shared<Matrix>(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, rdm3->data(), size, fock_->data(), 1, 0.0, work2->data(), 1);
  auto work = make_shared<Matrix>(dim, dim);
  sort_indices<0,2,1,3,0,1,1,1>(work2->data(), work->data(), nact, nact, nact, nact);
  work_xx_->copy_block(dim*jst, dim*ist, dim, dim, work);
}


void Denom::init_hh_(const int jst, const int ist, shared_ptr<const RDM<1>> rdm1, shared_ptr<const RDM<2>> rdm2, shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact;
  const size_t size = dim*dim;
  shared_ptr<RDM<2>> ovl = rdm2->copy();
  for (int i2 = 0; i2 != nact; ++i2)
    for (int i1 = 0; i1 != nact; ++i1)
      for (int i3 = 0; i3 != nact; ++i3)
        for (int i0 = 0; i0 != nact; ++i0) {
          double a = 0.0;
          if (i1 == i3)             a +=        rdm1->element(i0, i2);
          if (i1 == i2)             a += -2.0 * rdm1->element(i0, i3);
          if (i0 == i3)             a += -2.0 * rdm1->element(i1, i2);
          if (i0 == i2)             a +=        rdm1->element(i1, i3);
          if (i0 == i3 && i1 == i2) a +=  4.0;
          if (i0 == i2 && i1 == i3) a += -2.0;
          ovl->element(i0, i3, i1, i2) += a;
        }

  auto shalf = make_shared<Matrix>(dim, dim);
  sort_indices<0,2,1,3,0,1,1,1>(ovl->data(), shalf->data(), nact, nact, nact, nact);
  shalf_hh_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  shared_ptr<RDM<3>> r3 = rdm3->copy();
  for (int i2 = 0; i2 != nact; ++i2)
    for (int i3 = 0; i3 != nact; ++i3)
      for (int i4 = 0; i4 != nact; ++i4)
        for (int i1 = 0; i1 != nact; ++i1)
          for (int i5 = 0; i5 != nact; ++i5)
            for (int i0 = 0; i0 != nact; ++i0) {
              double a = 0.0;
              if (i3 == i5)                         a +=        rdm2->element(i1, i4, i0, i2);
              if (i3 == i4)                         a +=        rdm2->element(i0, i5, i1, i2);
              if (i1 == i5)                         a +=        rdm2->element(i0, i4, i3, i2);
              if (i1 == i5 && i3 == i4)             a +=        rdm1->element(i0, i2);
              if (i1 == i4)                         a += -2.0 * rdm2->element(i0, i5, i3, i2);
              if (i1 == i4 && i3 == i5)             a += -2.0 * rdm1->element(i0, i2);
              if (i1 == i2)                         a +=        rdm2->element(i0, i5, i3, i4);
              if (i1 == i2 && i3 == i5)             a +=        rdm1->element(i0, i4);
              if (i1 == i2 && i3 == i4)             a += -2.0 * rdm1->element(i0, i5);
              if (i0 == i5)                         a += -2.0 * rdm2->element(i1, i4, i3, i2);
              if (i0 == i5 && i3 == i4)             a += -2.0 * rdm1->element(i1, i2);
              if (i1 == i2 && i0 == i5)             a += -2.0 * rdm1->element(i3, i4);
              if (i1 == i2 && i0 == i5 && i3 == i4) a +=  4.0;
              if (i0 == i4)                         a +=        rdm2->element(i1, i5, i3, i2);
              if (i0 == i4 && i3 == i5)             a +=        rdm1->element(i1, i2);
              if (i1 == i2 && i0 == i4)             a +=        rdm1->element(i3, i5);
              if (i1 == i2 && i0 == i4 && i3 == i5) a += -2.0;
              if (i0 == i2)                         a +=        rdm2->element(i3, i5, i1, i4);
              if (i0 == i2 && i3 == i5)             a += -2.0 * rdm1->element(i1, i4);
              if (i0 == i2 && i3 == i4)             a +=        rdm1->element(i1, i5);
              if (i1 == i5 && i0 == i2)             a +=        rdm1->element(i3, i4);
              if (i0 == i2 && i1 == i5 && i3 == i4) a += -2.0;
              if (i0 == i2 && i1 == i4)             a += -2.0 * rdm1->element(i3, i5);
              if (i1 == i4 && i3 == i5 && i0 == i2) a +=  4.0;
              if (i0 == i5 && i1 == i4)             a +=  4.0 * rdm1->element(i3, i2);
              if (i0 == i4 && i1 == i5)             a += -2.0 * rdm1->element(i3, i2);
              r3->element(i0, i5, i1, i4, i3, i2) += a;
            }

  auto work2 = make_shared<Matrix>(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, r3->data(), size, fock_->data(), 1, 0.0, work2->data(), 1);
  auto work = make_shared<Matrix>(dim, dim);
  sort_indices<0,2,1,3,0,1,1,1>(work2->data(), work->data(), nact, nact, nact, nact);
  work_hh_->copy_block(dim*jst, dim*ist, dim, dim, work);
}



void Denom::init_xh_(const int jst, const int ist, shared_ptr<const RDM<1>> rdm1, shared_ptr<const RDM<2>> rdm2, shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact;
  const size_t size = dim*dim;
  auto shalf = make_shared<Matrix>(dim*2, dim*2);
  RDM<2> ovl1 = *rdm2;
  RDM<2> ovl4 = *rdm2;
  ovl1.scale(-1.0);
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k) {
        ovl1.element(i, i, k, j) += 2.0 * rdm1->element(k, j);
        ovl4.element(k, i, i, j) +=       rdm1->element(k, j);
      }

  Matrix work(dim, dim);
  sort_indices<3,0,2,1,0,1,1,1>(ovl1.data(), work.data(), nact, nact, nact, nact);
  shalf->add_block(1.0, dim, dim, dim, dim, work);

  sort_indices<0,1,3,2,0,1,2,1>(ovl4.data(), work.data(), nact, nact, nact, nact);
  shalf->add_block(1.0, 0, 0, dim, dim, work);

  shalf->add_block(-0.5, dim, 0, dim, dim, work);
  shalf->add_block(-0.5, 0, dim, dim, dim, work);

  shalf_xh_->copy_block(dim*jst*2, dim*ist*2, dim*2, dim*2, shalf);

  shared_ptr<RDM<3>> d0 = rdm3->copy();
  shared_ptr<RDM<3>> d3 = rdm3->copy();
  d0->scale(-1.0);
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i4 = 0; i4 != nact; ++i4)
      for (int i3 = 0; i3 != nact; ++i3)
        for (int i2 = 0; i2 != nact; ++i2)
          for (int i1 = 0; i1 != nact; ++i1)
            for (int i0 = 0; i0 != nact; ++i0) {
              {
                double a = 0.0;
                if (i3 == i4)             a += rdm2->element(i0, i1, i2, i5);
                if (i1 == i2)             a += rdm2->element(i0, i3, i4, i5);
                if (i1 == i2 && i3 == i4) a += rdm1->element(i0, i5);
                if (i1 == i4)             a += rdm2->element(i2, i3, i0, i5);
                d3->element(i0, i1, i4, i5, i2, i3) += a;
              } {
                double b = 0.0;
                if (i3 == i4)             b += -1.0 * rdm2->element(i2, i1, i0, i5);
                if (i1 == i2)             b += -1.0 * rdm2->element(i4, i3, i0, i5);
                if (i3 == i4 && i1 == i2) b +=  2.0 * rdm1->element(i0, i5);
                if (i1 == i4)             b +=  2.0 * rdm2->element(i2, i3, i0, i5);
                d0->element(i4, i1, i0, i5, i2, i3) += b;
              }
            }

  Matrix work2(dim, dim);

  dgemv_("N", size, nact*nact, 1.0, d0->data(), size, fock_->data(), 1, 0.0, work2.data(), 1);
  sort_indices<3,0,2,1,0,1,1,1>(work2.data(), work.data(), nact, nact, nact, nact);
  auto num = make_shared<Matrix>(dim*2, dim*2);
  num->add_block(1.0, dim, dim, dim, dim, work);

  dgemv_("N", size, nact*nact, 1.0, d3->data(), size, fock_->data(), 1, 0.0, work2.data(), 1);
  sort_indices<0,1,3,2,0,1,2,1>(work2.data(), work.data(), nact, nact, nact, nact);
  num->add_block(1.0, 0, 0, dim, dim, work);
  num->add_block(-0.5, dim, 0, dim, dim, work);
  num->add_block(-0.5, 0, dim, dim, dim, work);

  work_xh_->copy_block(dim*jst*2, dim*ist*2, dim*2, dim*2, num);
}


void Denom::init_xhh_(const int jst, const int ist, shared_ptr<const RDM<1>> rdm1, shared_ptr<const RDM<2>> rdm2, shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact*nact;
  const size_t size = dim*dim;
  shared_ptr<RDM<3>> ovl = rdm3->copy();
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i0 = 0; i0 != nact; ++i0)
      for (int i4 = 0; i4 != nact; ++i4)
        for (int i3 = 0; i3 != nact; ++i3)
          for (int i2 = 0; i2 != nact; ++i2)
            for (int i1 = 0; i1 != nact; ++i1) {
              if (i2 == i3) ovl->element(i1, i2, i3, i4, i0, i5) += 1.0 * rdm2->element(i1, i4, i0, i5);
            }
  auto shalf = make_shared<Matrix>(dim, dim);
  sort_indices<4,0,1,5,3,2,0,1,1,1>(ovl->data(), shalf->data(), nact, nact, nact, nact, nact, nact);
  shalf_xhh_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  shared_ptr<RDM<4>> r4 = rdm4->copy();
  for (int i4 = 0; i4 != nact; ++i4)
    for (int i3 = 0; i3 != nact; ++i3)
      for (int i7 = 0; i7 != nact; ++i7)
        for (int i0 = 0; i0 != nact; ++i0)
          for (int i6 = 0; i6 != nact; ++i6)
            for (int i5 = 0; i5 != nact; ++i5)
              for (int i2 = 0; i2 != nact; ++i2)
                for (int i1 = 0; i1 != nact; ++i1) {
                  double a = 0.0;
                  if (i4 == i5)             a += 1.0 * rdm3->element(i1, i2, i3, i6, i0, i7);
                  if (i2 == i3)             a += 1.0 * rdm3->element(i1, i4, i5, i6, i0, i7);
                  if (i2 == i3 && i4 == i5) a += 1.0 * rdm2->element(i1, i6, i0, i7);
                  if (i2 == i5)             a += 1.0 * rdm3->element(i3, i4, i1, i6, i0, i7);
                  r4->element(i1, i2, i5, i6, i0, i7, i3, i4) += a;
                }
  Matrix work2(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, r4->data(), size, fock_->data(), 1, 0.0, work2.data(), 1);
  auto fss = make_shared<Matrix>(dim, dim);
  sort_indices<4,0,1,5,3,2,0,1,1,1>(work2.data(), fss->data(), nact, nact, nact, nact, nact, nact);
  work_xhh_->copy_block(dim*jst, dim*ist, dim, dim, fss);
}


void Denom::init_xxh_(const int jst, const int ist, shared_ptr<const RDM<1>> rdm1, shared_ptr<const RDM<2>> rdm2, shared_ptr<const RDM<3>> rdm3, shared_ptr<const RDM<4>> rdm4) {
  const size_t nact = rdm1->norb();
  const size_t dim  = nact*nact*nact;
  const size_t size = dim*dim;
  shared_ptr<RDM<3>> ovl = rdm3->copy();
  ovl->scale(-1.0);
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i4 = 0; i4 != nact; ++i4)
      for (int i2 = 0; i2 != nact; ++i2)
        for (int i3 = 0; i3 != nact; ++i3)
          for (int i1 = 0; i1 != nact; ++i1)
            for (int i0 = 0; i0 != nact; ++i0) {
              double a = 0.0;
              if (i2 == i4)             a += -1.0 * rdm2->element(i0, i1, i3, i5);
              if (i2 == i3)             a +=  2.0 * rdm2->element(i0, i1, i4, i5);
              if (i1 == i4)             a += -1.0 * rdm2->element(i3, i2, i0, i5);
              if (i1 == i4 && i2 == i3) a +=  2.0 * rdm1->element(i0, i5);
              if (i1 == i3)             a += -1.0 * rdm2->element(i0, i2, i4, i5);
              if (i1 == i3 && i2 == i4) a += -1.0 * rdm1->element(i0, i5);
              ovl->element(i0, i1, i3, i2, i4, i5) += a;
            }
  auto shalf = make_shared<Matrix>(dim, dim);
  sort_indices<0,1,3,5,4,2,0,1,1,1>(ovl->data(), shalf->data(), nact, nact, nact, nact, nact, nact);
  shalf_xxh_->copy_block(dim*jst, dim*ist, dim, dim, shalf);

  shared_ptr<RDM<4>> r4 = rdm4->copy();
  r4->scale(-1.0);
  for (int i4 = 0; i4 != nact; ++i4)
    for (int i3 = 0; i3 != nact; ++i3)
      for (int i7 = 0; i7 != nact; ++i7)
        for (int i6 = 0; i6 != nact; ++i6)
          for (int i2 = 0; i2 != nact; ++i2)
            for (int i5 = 0; i5 != nact; ++i5)
              for (int i1 = 0; i1 != nact; ++i1)
                for (int i0 = 0; i0 != nact; ++i0) {
                  double a = 0.0;
                  if (i4 == i6)                         a += -1.0 * rdm3->element(i0, i1, i5, i2, i3, i7);
                  if (i4 == i5)                         a += -1.0 * rdm3->element(i0, i1, i3, i2, i6, i7);
                  if (i2 == i6)                         a += -1.0 * rdm3->element(i0, i1, i3, i4, i5, i7);
                  if (i2 == i6 && i4 == i5)             a += -1.0 * rdm2->element(i0, i1, i3, i7);
                  if (i2 == i5)                         a +=  2.0 * rdm3->element(i0, i1, i3, i4, i6, i7);
                  if (i2 == i5 && i4 == i6)             a +=  2.0 * rdm2->element(i0, i1, i3, i7);
                  if (i2 == i3)                         a += -1.0 * rdm3->element(i0, i1, i5, i4, i6, i7);
                  if (i2 == i3 && i4 == i6)             a += -1.0 * rdm2->element(i0, i1, i5, i7);
                  if (i2 == i3 && i4 == i5)             a +=  2.0 * rdm2->element(i0, i1, i6, i7);
                  if (i1 == i6)                         a += -1.0 * rdm3->element(i3, i4, i5, i2, i0, i7);
                  if (i1 == i6 && i4 == i5)             a += -1.0 * rdm2->element(i3, i2, i0, i7);
                  if (i1 == i6 && i2 == i3)             a += -1.0 * rdm2->element(i5, i4, i0, i7);
                  if (i1 == i6 && i2 == i3 && i4 == i5) a +=  2.0 * rdm1->element(i0, i7);
                  if (i1 == i5)                         a += -1.0 * rdm3->element(i0, i2, i3, i4, i6, i7);
                  if (i1 == i5 && i4 == i6)             a += -1.0 * rdm2->element(i0, i2, i3, i7);
                  if (i2 == i3 && i1 == i5)             a += -1.0 * rdm2->element(i0, i4, i6, i7);
                  if (i2 == i3 && i1 == i5 && i4 == i6) a += -1.0 * rdm1->element(i0, i7);
                  if (i1 == i3)                         a += -1.0 * rdm3->element(i0, i4, i5, i2, i6, i7);
                  if (i1 == i3 && i4 == i6)             a += -1.0 * rdm2->element(i5, i2, i0, i7);
                  if (i1 == i3 && i4 == i5)             a += -1.0 * rdm2->element(i0, i2, i6, i7);
                  if (i2 == i6 && i1 == i3)             a += -1.0 * rdm2->element(i0, i4, i5, i7);
                  if (i2 == i6 && i1 == i3 && i4 == i5) a += -1.0 * rdm1->element(i0, i7);
                  if (i1 == i3 && i2 == i5)             a +=  2.0 * rdm2->element(i0, i4, i6, i7);
                  if (i4 == i6 && i2 == i5 && i1 == i3) a +=  2.0 * rdm1->element(i0, i7);
                  if (i1 == i6 && i2 == i5)             a +=  2.0 * rdm2->element(i3, i4, i0, i7);
                  if (i1 == i5 && i2 == i6)             a += -1.0 * rdm2->element(i3, i4, i0, i7);
                  r4->element(i0, i1, i5, i2, i6, i7, i3, i4) += a;
                }
  Matrix work2(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, r4->data(), size, fock_->data(), 1, 0.0, work2.data(), 1);
  auto fss = make_shared<Matrix>(dim, dim);
  sort_indices<0,1,3,5,4,2,0,1,1,1>(work2.data(), fss->data(), nact, nact, nact, nact, nact, nact);
  work_xxh_->copy_block(dim*jst, dim*ist, dim, dim, fss);
}

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


#include <src/smith/denom.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;
using namespace SMITH;

Denom::Denom(const RDM<1>& rdm1, const RDM<2>& rdm2, const RDM<3>& rdm3, const RDM<4>& rdm4, const Matrix& fock) {
  init_hh_(rdm1, rdm2, rdm3, rdm4, fock);
  init_xh_(rdm1, rdm2, rdm3, rdm4, fock);
  init_xhh_(rdm1, rdm2, rdm3, rdm4, fock);
  init_xxh_(rdm1, rdm2, rdm3, rdm4, fock);
}


void Denom::init_hh_(const RDM<1>& rdm1, const RDM<2>& rdm2, const RDM<3>& rdm3, const RDM<4>& rdm4, const Matrix& fock) {
  const size_t nact = rdm1.norb();
  const size_t dim  = nact*nact;
  const size_t size = dim*dim;
  RDM<2> ovl = rdm2;
  for (int i2 = 0; i2 != nact; ++i2)
    for (int i1 = 0; i1 != nact; ++i1)
      for (int i3 = 0; i3 != nact; ++i3)
        for (int i0 = 0; i0 != nact; ++i0) {
          if (i1 == i3)             ovl.element(i0, i3, i1, i2) +=        rdm1.element(i0, i2);
          if (i1 == i2)             ovl.element(i0, i3, i1, i2) += -2.0 * rdm1.element(i0, i3);
          if (i0 == i3)             ovl.element(i0, i3, i1, i2) += -2.0 * rdm1.element(i1, i2);
          if (i0 == i2)             ovl.element(i0, i3, i1, i2) +=        rdm1.element(i1, i3);
          if (i0 == i3 && i1 == i2) ovl.element(i0, i3, i1, i2) +=  4.0;
          if (i0 == i2 && i1 == i3) ovl.element(i0, i3, i1, i2) += -2.0;
        }

  Matrix shalf(dim, dim);
  sort_indices<0,2,1,3,0,1,1,1>(ovl.data(), shalf.data(), nact, nact, nact, nact);
  shalf.inverse_half(1.0e-9);

  RDM<3> r3 = rdm3;
  for (int i2 = 0; i2 != nact; ++i2)
    for (int i3 = 0; i3 != nact; ++i3)
      for (int i4 = 0; i4 != nact; ++i4)
        for (int i1 = 0; i1 != nact; ++i1)
          for (int i5 = 0; i5 != nact; ++i5)
            for (int i0 = 0; i0 != nact; ++i0) {
              if (i3 == i5)                         r3.element(i0, i5, i1, i4, i3, i2) +=        rdm2.element(i1, i4, i0, i2);
              if (i3 == i4)                         r3.element(i0, i5, i1, i4, i3, i2) +=        rdm2.element(i0, i5, i1, i2);
              if (i1 == i5)                         r3.element(i0, i5, i1, i4, i3, i2) +=        rdm2.element(i0, i4, i3, i2);
              if (i1 == i5 && i3 == i4)             r3.element(i0, i5, i1, i4, i3, i2) +=        rdm1.element(i0, i2);
              if (i1 == i4)                         r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm2.element(i0, i5, i3, i2);
              if (i1 == i4 && i3 == i5)             r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm1.element(i0, i2);
              if (i1 == i2)                         r3.element(i0, i5, i1, i4, i3, i2) +=        rdm2.element(i0, i5, i3, i4);
              if (i1 == i2 && i3 == i5)             r3.element(i0, i5, i1, i4, i3, i2) +=        rdm1.element(i0, i4);
              if (i1 == i2 && i3 == i4)             r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm1.element(i0, i5);
              if (i0 == i5)                         r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm2.element(i1, i4, i3, i2);
              if (i0 == i5 && i3 == i4)             r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm1.element(i1, i2);
              if (i1 == i2 && i0 == i5)             r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm1.element(i3, i4);
              if (i1 == i2 && i0 == i5 && i3 == i4) r3.element(i0, i5, i1, i4, i3, i2) +=  4.0;
              if (i0 == i4)                         r3.element(i0, i5, i1, i4, i3, i2) +=        rdm2.element(i1, i5, i3, i2);
              if (i0 == i4 && i3 == i5)             r3.element(i0, i5, i1, i4, i3, i2) +=        rdm1.element(i1, i2);
              if (i1 == i2 && i0 == i4)             r3.element(i0, i5, i1, i4, i3, i2) +=        rdm1.element(i3, i5);
              if (i1 == i2 && i0 == i4 && i3 == i5) r3.element(i0, i5, i1, i4, i3, i2) += -2.0;
              if (i0 == i2)                         r3.element(i0, i5, i1, i4, i3, i2) +=        rdm2.element(i3, i5, i1, i4);
              if (i0 == i2 && i3 == i5)             r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm1.element(i1, i4);
              if (i0 == i2 && i3 == i4)             r3.element(i0, i5, i1, i4, i3, i2) +=        rdm1.element(i1, i5);
              if (i1 == i5 && i0 == i2)             r3.element(i0, i5, i1, i4, i3, i2) +=        rdm1.element(i3, i4);
              if (i0 == i2 && i1 == i5 && i3 == i4) r3.element(i0, i5, i1, i4, i3, i2) += -2.0;
              if (i0 == i2 && i1 == i4)             r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm1.element(i3, i5);
              if (i1 == i4 && i3 == i5 && i0 == i2) r3.element(i0, i5, i1, i4, i3, i2) +=  4.0;
              if (i0 == i5 && i1 == i4)             r3.element(i0, i5, i1, i4, i3, i2) +=  4.0 * rdm1.element(i3, i2);
              if (i0 == i4 && i1 == i5)             r3.element(i0, i5, i1, i4, i3, i2) += -2.0 * rdm1.element(i3, i2);
            }

  Matrix work(dim, dim);
  Matrix work2(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, r3.data(), size, fock.data(), 1, 0.0, work2.data(), 1);

  sort_indices<0,2,1,3,0,1,1,1>(work2.data(), work.data(), nact, nact, nact, nact);
  Matrix fss = shalf % work * shalf;
  denom_hh_ = unique_ptr<double[]>(new double[dim]);
  fss.diagonalize(denom_hh_.get());
  shalf_hh_ = shared_ptr<const Matrix>(new Matrix(fss % shalf));
}



void Denom::init_xh_(const RDM<1>& rdm1, const RDM<2>& rdm2, const RDM<3>& rdm3, const RDM<4>& rdm4, const Matrix& fock) {
  const size_t nact = rdm1.norb();
  const size_t dim  = nact*nact;
  const size_t size = dim*dim;
  Matrix shalf(dim*2, dim*2);
  RDM<2> ovl1 = rdm2;
  RDM<2> ovl4 = rdm2;
  ovl1.scale(-1.0);
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k) {
        ovl1.element(i, i, k, j) += 2.0 * rdm1.element(k, j);
        ovl4.element(k, i, i, j) +=       rdm1.element(k, j);
      }

  Matrix work(dim, dim);
  sort_indices<3,0,2,1,0,1,1,1>(ovl1.data(), work.data(), nact, nact, nact, nact);
  shalf.add_block(dim, dim, dim, dim, work);

  sort_indices<0,1,3,2,0,1,2,1>(ovl4.data(), work.data(), nact, nact, nact, nact);
  shalf.add_block(0, 0, dim, dim, work);

  work.scale(-0.5);
  shalf.add_block(dim, 0, dim, dim, work);
  shalf.add_block(0, dim, dim, dim, work);

  shalf.inverse_half(1.0e-9);

  RDM<3> d0 = rdm3;
  RDM<3> d3 = rdm3; 
  d0.scale(-1.0);
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i4 = 0; i4 != nact; ++i4)
      for (int i3 = 0; i3 != nact; ++i3)
        for (int i2 = 0; i2 != nact; ++i2)
          for (int i1 = 0; i1 != nact; ++i1)
            for (int i0 = 0; i0 != nact; ++i0) {
              if (i3 == i4)             d3.element(i0, i1, i4, i5, i2, i3) += rdm2.element(i0, i1, i2, i5);
              if (i1 == i2)             d3.element(i0, i1, i4, i5, i2, i3) += rdm2.element(i0, i3, i4, i5);
              if (i1 == i2 && i3 == i4) d3.element(i0, i1, i4, i5, i2, i3) += rdm1.element(i0, i5);
              if (i1 == i4)             d3.element(i0, i1, i4, i5, i2, i3) += rdm2.element(i2, i3, i0, i5);
              if (i3 == i4)             d0.element(i4, i1, i0, i5, i2, i3) += -1.0 * rdm2.element(i2, i1, i0, i5);
              if (i1 == i2)             d0.element(i4, i1, i0, i5, i2, i3) += -1.0 * rdm2.element(i4, i3, i0, i5);
              if (i3 == i4 && i1 == i2) d0.element(i4, i1, i0, i5, i2, i3) +=  2.0 * rdm1.element(i0, i5);
              if (i1 == i4)             d0.element(i4, i1, i0, i5, i2, i3) +=  2.0 * rdm2.element(i2, i3, i0, i5);
            }

  Matrix work2(dim, dim);

  dgemv_("N", size, nact*nact, 1.0, d0.data(), size, fock.data(), 1, 0.0, work2.data(), 1);
  sort_indices<3,0,2,1,0,1,1,1>(work2.data(), work.data(), nact, nact, nact, nact);
  Matrix num(dim*2, dim*2);
  num.add_block(dim, dim, dim, dim, work);

  dgemv_("N", size, nact*nact, 1.0, d3.data(), size, fock.data(), 1, 0.0, work2.data(), 1);
  sort_indices<0,1,3,2,0,1,2,1>(work2.data(), work.data(), nact, nact, nact, nact);
  num.add_block(0, 0, dim, dim, work);
  work.scale(-0.5);
  num.add_block(dim, 0, dim, dim, work);
  num.add_block(0, dim, dim, dim, work);

  Matrix fss = shalf % num * shalf;
  denom_xh_ = unique_ptr<double[]>(new double[2*dim]);
  fss.diagonalize(denom_xh_.get());
  shalf_xh_ = shared_ptr<const Matrix>(new Matrix(fss % shalf));
}


void Denom::init_xhh_(const RDM<1>& rdm1, const RDM<2>& rdm2, const RDM<3>& rdm3, const RDM<4>& rdm4, const Matrix& fock) {
  const size_t nact = rdm1.norb();
  const size_t dim  = nact*nact*nact;
  const size_t size = dim*dim;
  RDM<3> ovl = rdm3;
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i0 = 0; i0 != nact; ++i0)
      for (int i4 = 0; i4 != nact; ++i4)
        for (int i3 = 0; i3 != nact; ++i3)
          for (int i2 = 0; i2 != nact; ++i2)
            for (int i1 = 0; i1 != nact; ++i1) {
              if (i2 == i3) ovl.element(i1, i2, i3, i4, i0, i5) += 1.0 * rdm2.element(i1, i4, i0, i5);
            }
  Matrix shalf(dim, dim);
  sort_indices<4,0,1,5,3,2,0,1,1,1>(ovl.data(), shalf.data(), nact, nact, nact, nact, nact, nact);
  shalf.inverse_half(1.0e-9);

  RDM<4> r4 = rdm4;
  for (int i4 = 0; i4 != nact; ++i4)
    for (int i3 = 0; i3 != nact; ++i3)
      for (int i7 = 0; i7 != nact; ++i7)
        for (int i0 = 0; i0 != nact; ++i0)
          for (int i6 = 0; i6 != nact; ++i6)
            for (int i5 = 0; i5 != nact; ++i5)
              for (int i2 = 0; i2 != nact; ++i2)
                for (int i1 = 0; i1 != nact; ++i1) {
                  if (i4 == i5)             r4.element(i1, i2, i5, i6, i0, i7, i3, i4) += 1.0 * rdm3.element(i1, i2, i3, i6, i0, i7);
                  if (i2 == i3)             r4.element(i1, i2, i5, i6, i0, i7, i3, i4) += 1.0 * rdm3.element(i1, i4, i5, i6, i0, i7);
                  if (i2 == i3 && i4 == i5) r4.element(i1, i2, i5, i6, i0, i7, i3, i4) += 1.0 * rdm2.element(i1, i6, i0, i7);
                  if (i2 == i5)             r4.element(i1, i2, i5, i6, i0, i7, i3, i4) += 1.0 * rdm3.element(i3, i4, i1, i6, i0, i7);
                }
  Matrix work2(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, r4.data(), size, fock.data(), 1, 0.0, work2.data(), 1);
  Matrix fss(dim, dim);
  sort_indices<4,0,1,5,3,2,0,1,1,1>(work2.data(), fss.data(), nact, nact, nact, nact, nact, nact);
  fss = shalf % fss * shalf;
  denom_xhh_ = unique_ptr<double[]>(new double[dim]);
  fss.diagonalize(denom_xhh_.get());
  shalf_xhh_ = shared_ptr<const Matrix>(new Matrix(fss % shalf));
}


void Denom::init_xxh_(const RDM<1>& rdm1, const RDM<2>& rdm2, const RDM<3>& rdm3, const RDM<4>& rdm4, const Matrix& fock) {
  const size_t nact = rdm1.norb();
  const size_t dim  = nact*nact*nact;
  const size_t size = dim*dim;
  RDM<3> ovl = rdm3;
  ovl.scale(-1.0);
  for (int i5 = 0; i5 != nact; ++i5)
    for (int i4 = 0; i4 != nact; ++i4)
      for (int i2 = 0; i2 != nact; ++i2)
        for (int i3 = 0; i3 != nact; ++i3)
          for (int i1 = 0; i1 != nact; ++i1)
            for (int i0 = 0; i0 != nact; ++i0) {
              if (i2 == i4)             ovl.element(i0, i1, i3, i2, i4, i5) += -1.0 * rdm2.element(i0, i1, i3, i5);
              if (i2 == i3)             ovl.element(i0, i1, i3, i2, i4, i5) +=  2.0 * rdm2.element(i0, i1, i4, i5);
              if (i1 == i4)             ovl.element(i0, i1, i3, i2, i4, i5) += -1.0 * rdm2.element(i3, i2, i0, i5);
              if (i1 == i4 && i2 == i3) ovl.element(i0, i1, i3, i2, i4, i5) +=  2.0 * rdm1.element(i0, i5);
              if (i1 == i3)             ovl.element(i0, i1, i3, i2, i4, i5) += -1.0 * rdm2.element(i0, i2, i4, i5);
              if (i1 == i3 && i2 == i4) ovl.element(i0, i1, i3, i2, i4, i5) += -1.0 * rdm1.element(i0, i5);
            }
  Matrix shalf(dim, dim);
  sort_indices<0,1,3,5,4,2,0,1,1,1>(ovl.data(), shalf.data(), nact, nact, nact, nact, nact, nact);
  shalf.inverse_half(1.0e-9);

  RDM<4> r4 = rdm4;
  r4.scale(-1.0);
  for (int i4 = 0; i4 != nact; ++i4)
    for (int i3 = 0; i3 != nact; ++i3)
      for (int i7 = 0; i7 != nact; ++i7)
        for (int i6 = 0; i6 != nact; ++i6)
          for (int i2 = 0; i2 != nact; ++i2)
            for (int i5 = 0; i5 != nact; ++i5)
              for (int i1 = 0; i1 != nact; ++i1)
                for (int i0 = 0; i0 != nact; ++i0) {
                  if (i4 == i6)                         r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3.element(i0, i1, i5, i2, i3, i7);
                  if (i4 == i5)                         r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3.element(i0, i1, i3, i2, i6, i7);
                  if (i2 == i6)                         r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3.element(i0, i1, i3, i4, i5, i7);
                  if (i2 == i6 && i4 == i5)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i0, i1, i3, i7);
                  if (i2 == i5)                         r4.element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm3.element(i0, i1, i3, i4, i6, i7);
                  if (i2 == i5 && i4 == i6)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm2.element(i0, i1, i3, i7);
                  if (i2 == i3)                         r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3.element(i0, i1, i5, i4, i6, i7);
                  if (i2 == i3 && i4 == i6)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i0, i1, i5, i7);
                  if (i2 == i3 && i4 == i5)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm2.element(i0, i1, i6, i7);
                  if (i1 == i6)                         r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3.element(i3, i4, i5, i2, i0, i7);
                  if (i1 == i6 && i4 == i5)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i3, i2, i0, i7);
                  if (i1 == i6 && i2 == i3)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i5, i4, i0, i7);
                  if (i1 == i6 && i2 == i3 && i4 == i5) r4.element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm1.element(i0, i7);
                  if (i1 == i5)                         r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3.element(i0, i2, i3, i4, i6, i7);
                  if (i1 == i5 && i4 == i6)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i0, i2, i3, i7);
                  if (i2 == i3 && i1 == i5)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i0, i4, i6, i7);
                  if (i2 == i3 && i1 == i5 && i4 == i6) r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm1.element(i0, i7);
                  if (i1 == i3)                         r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm3.element(i0, i4, i5, i2, i6, i7);
                  if (i1 == i3 && i4 == i6)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i5, i2, i0, i7);
                  if (i1 == i3 && i4 == i5)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i0, i2, i6, i7);
                  if (i2 == i6 && i1 == i3)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i0, i4, i5, i7);
                  if (i2 == i6 && i1 == i3 && i4 == i5) r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm1.element(i0, i7);
                  if (i1 == i3 && i2 == i5)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm2.element(i0, i4, i6, i7);
                  if (i4 == i6 && i2 == i5 && i1 == i3) r4.element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm1.element(i0, i7);
                  if (i1 == i6 && i2 == i5)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) +=  2.0 * rdm2.element(i3, i4, i0, i7);
                  if (i1 == i5 && i2 == i6)             r4.element(i0, i1, i5, i2, i6, i7, i3, i4) += -1.0 * rdm2.element(i3, i4, i0, i7);
                }
  Matrix work2(dim, dim);
  dgemv_("N", size, nact*nact, 1.0, r4.data(), size, fock.data(), 1, 0.0, work2.data(), 1);
  Matrix fss(dim, dim);
  sort_indices<0,1,3,5,4,2,0,1,1,1>(work2.data(), fss.data(), nact, nact, nact, nact, nact, nact);
  fss = shalf % fss * shalf;
  denom_xxh_ = unique_ptr<double[]>(new double[dim]);
  fss.diagonalize(denom_xxh_.get());
  shalf_xxh_ = shared_ptr<const Matrix>(new Matrix(fss % shalf)); 
}


//
// BAGEL - Parallel electron correlation program.
// Filename: nevpt2_rdm.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/nevpt2/nevpt2.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;


tuple<shared_ptr<const Matrix>, shared_ptr<const Matrix>, shared_ptr<const Matrix>, shared_ptr<const Matrix>, shared_ptr<const Matrix>>
  NEVPT2::compute_asrdm(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> rdm2, shared_ptr<const Matrix> rdm3, shared_ptr<const Matrix> rdm4) const {

  const int nact = rdm1->ndim();

  static auto id2 = [&nact](                          const int k, const int l) { return        (       (k+nact*l)); };
  static auto id3 = [&nact](             const int j, const int k, const int l) { return        (j+nact*(k+nact*l)); };
  static auto id4 = [&nact](const int i, const int j, const int k, const int l) { return i+nact*(j+nact*(k+nact*l)); };

  shared_ptr<Matrix> srdm2 = rdm2->clone(); // S(a,b,c,d) = <0|a+p bp cq d+q|0>
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        for (int l = 0; l != nact; ++l)
          srdm2->element(l+nact*k,j+nact*i) = -rdm2->element(l+nact*i,k+nact*j) + (i == j ? 2.0*rdm1->element(l,k) : 0.0) - (i == k ? rdm1->element(l,j) : 0.0);
  // <a+ a b+ b> and <a+ a b+ b c+ c>
  shared_ptr<Matrix> ardm2 = rdm2->clone();
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k) {
        for (int l = 0; l != nact; ++l)
          ardm2->element(l+nact*k,j+nact*i) += rdm2->element(l+nact*j,k+nact*i);
        ardm2->element(k+nact*j,j+nact*i) += rdm1->element(k,i);
      }
  shared_ptr<Matrix> ardm3 = rdm3->clone();
  shared_ptr<Matrix> srdm3 = rdm3->clone(); // <a+ a b b+ c+ c>
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        for (int l = 0; l != nact; ++l)
          for (int m = 0; m != nact; ++m) {
            for (int n = 0; n != nact; ++n)
              ardm3->element(id3(n,m,l),id3(k,j,i)) += rdm3->element(id3(n,l,j),id3(m,k,i));
            ardm3->element(id3(m,l,l),id3(k,j,i)) += ardm2->element(m+nact*k,j+nact*i);
            ardm3->element(id3(m,l,k),id3(j,j,i)) += rdm2->element(m+nact*k,l+nact*i);
            ardm3->element(id3(m,l,k),id3(j,l,i)) += rdm2->element(m+nact*k,i+nact*j);

            srdm3->element(id3(m,l,k),id3(k,j,i)) += 2.0*ardm2->element(id2(m,l),id2(j,i));
          }
  SMITH::sort_indices<0,2,1,3,1,1,-1,1>(ardm3->data(), srdm3->data(), nact*nact, nact, nact, nact*nact);
  shared_ptr<Matrix> ardm4 = rdm4->clone();
  for (int h = 0; h != nact; ++h)
    for (int g = 0; g != nact; ++g)
      for (int f = 0; f != nact; ++f)
        for (int e = 0; e != nact; ++e)
          for (int d = 0; d != nact; ++d)
            for (int c = 0; c != nact; ++c)
              for (int b = 0; b != nact; ++b)
                for (int a = 0; a != nact; ++a) {
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == c ? 1.0 : 0.0) * ardm3->element(id3(a,d,e),id3(f,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) -= (d == e && b == c ? 1.0 : 0.0) * ardm2->element(id2(a,f),id2(g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (d == e ? 1.0 : 0.0) * ardm3->element(id3(a,b,c),id3(f,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) -= (b == e && c == f ? 1.0 : 0.0) * ardm2->element(id2(a,d),id2(g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == e ? 1.0 : 0.0) * ardm3->element(id3(a,f,c),id3(d,g,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (f == g ? 1.0 : 0.0) *  rdm3->element(id3(a,c,e),id3(b,d,h));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (d == g ? 1.0 : 0.0) *  rdm3->element(id3(a,c,e),id3(b,h,f));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += (b == g ? 1.0 : 0.0) *  rdm3->element(id3(a,c,e),id3(h,d,f));
                  ardm4->element(id4(a,b,c,d),id4(e,f,g,h)) += rdm4->element(id4(a,c,e,g),id4(b,d,f,h));
                }
  return make_tuple(ardm2, ardm3, ardm4, srdm2, srdm3);
}


tuple<shared_ptr<const Matrix>, shared_ptr<const Matrix>, shared_ptr<const Matrix>>
  NEVPT2::compute_hrdm(shared_ptr<const Matrix> rdm1, shared_ptr<const Matrix> rdm2, shared_ptr<const Matrix> rdm3, shared_ptr<const Matrix> srdm2) const {

  const int nact = rdm1->ndim();
  static auto id3 = [&nact](             const int j, const int k, const int l) { return        (j+nact*(k+nact*l)); };

  shared_ptr<Matrix> unit = rdm1->clone(); unit->unit();
  shared_ptr<const Matrix> hrdm1 = make_shared<Matrix>(*unit*2.0 - *rdm1);
  shared_ptr<Matrix> hrdm2 = rdm2->copy();
  for (int i = 0; i != nact; ++i) {
    for (int j = 0; j != nact; ++j) {
      for (int k = 0; k != nact; ++k) {
        hrdm2->element(j+nact*k, i+nact*k) += 2.0 * hrdm1->element(j,i);
        hrdm2->element(k+nact*j, i+nact*k) -= hrdm1->element(j,i);
        hrdm2->element(j+nact*k, k+nact*i) +=  rdm1->element(i,j);
        hrdm2->element(k+nact*j, k+nact*i) -= 2.0 *  rdm1->element(i,j);
      }
    }
  }
  shared_ptr<Matrix> hrdm3 = make_shared<Matrix>(*rdm3 * (-1.0));
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        for (int l = 0; l != nact; ++l)
          for (int m = 0; m != nact; ++m) {
            hrdm3->element(id3(l,k,m),id3(j,i,m)) += 2.0*hrdm2->element(l+nact*k,j+nact*i);
            hrdm3->element(id3(l,m,k),id3(j,i,m)) -=     hrdm2->element(l+nact*k,j+nact*i);
            hrdm3->element(id3(m,l,k),id3(j,i,m)) -=     hrdm2->element(l+nact*k,i+nact*j);
            hrdm3->element(id3(l,k,m),id3(j,m,i)) +=     srdm2->element(i+nact*k,l+nact*j);
            hrdm3->element(id3(l,m,k),id3(j,m,i)) -= 2.0*srdm2->element(i+nact*k,l+nact*j);
            hrdm3->element(id3(m,l,k),id3(j,m,i)) +=     srdm2->element(i+nact*k,l+nact*j);
            hrdm3->element(id3(l,k,m),id3(m,j,i)) -=      rdm2->element(i+nact*j,l+nact*k);
            hrdm3->element(id3(l,m,k),id3(m,j,i)) -=      rdm2->element(i+nact*j,k+nact*l);
            hrdm3->element(id3(m,l,k),id3(m,j,i)) += 2.0* rdm2->element(i+nact*j,k+nact*l);
          }
  return make_tuple(hrdm1, hrdm2, hrdm3);
}

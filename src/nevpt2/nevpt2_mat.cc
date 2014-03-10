//
// BAGEL - Parallel electron correlation program.
// Filename: nevpt2_mat.cc
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
#include <src/casscf/qvec.h>

using namespace std;
using namespace bagel;

void NEVPT2::compute_ints() {
  shared_ptr<const DFFullDist> full = casscf_->fci()->jop()->mo2e_1ext()->compute_second_transform(acoeff_)->apply_J();
  // integrals (ij|kl) and <ik|jl>
  shared_ptr<const Matrix> ints = full->form_4index(full, 1.0);
  auto tmp = make_shared<Matrix>(nact_*nact_, nact_*nact_);
  SMITH::sort_indices<0,2,1,3,0,1,1,1>(ints->data(), tmp->data(), nact_, nact_, nact_, nact_);
  ints2_ = tmp;
}


void NEVPT2::compute_kmat() {
  {
    auto kmat = make_shared<Matrix>(*fockact_c_ * *rdm1_);
    *kmat += Qvec(nact_, nact_, acoeff_, /*nclosed_*/0, casscf_->fci(), casscf_->fci()->rdm2(istate_));

    auto kmatp = make_shared<Matrix>(*kmat * (-1.0));
    *kmatp += *fockact_ * 2.0;

    kmat_ = kmat;
    kmatp_ = kmatp;
  }
  {
    auto compute_kmat = [](const int nact_, shared_ptr<const Matrix> rdm2, shared_ptr<const Matrix> rdm3, shared_ptr<const Matrix> fock,
                           shared_ptr<const Matrix> ints, const double sign) {
      auto out = rdm2->clone();
      // temp area
      shared_ptr<Matrix> four  = rdm2->clone();
      shared_ptr<Matrix> four2 = rdm2->clone();
      shared_ptr<Matrix> six   = rdm3->clone();

      // + (h)rdm2(i,j,k,m) * fockact_h_(m,l)
      dgemm_("N", "N", nact_*nact_*nact_, nact_, nact_, 1.0, rdm2->data(), nact_*nact_*nact_, fock->data(), nact_, 1.0, out->data(), nact_*nact_*nact_);
      // + (h)rdm2(i,j,m,k) * fockact_h_(m,l)
      SMITH::sort_indices<0,1,3,2,0,1,1,1>(rdm2->data(), four->data(), nact_, nact_, nact_, nact_);
      dgemm_("N", "N", nact_*nact_*nact_, nact_, nact_, 1.0, four->data(), nact_*nact_*nact_, fock->data(), nact_, 0.0, four2->data(), nact_*nact_*nact_);
      SMITH::sort_indices<0,1,3,2,1,1,1,1>(four2->data(), out->data(), nact_, nact_, nact_, nact_);
      // +/- (h)rdm2(i,j,m,n) * <mn|kl>
      dgemm_("N", "N", nact_*nact_, nact_*nact_, nact_*nact_, sign, rdm2->data(), nact_*nact_, ints->data(), nact_*nact_, 1.0, out->data(), nact_*nact_);
      // +/- (h)rdm3(i,j,x,y,l,z) * (<kx|yz>*)
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(rdm3->data(), six->data(), nact_*nact_, nact_*nact_, nact_, nact_); // (i,j,l,x,y,w)
      dgemm_("N", "T", nact_*nact_*nact_, nact_, nact_*nact_*nact_, sign, six->data(), nact_*nact_*nact_, ints->data(), nact_, 0.0, four->data(), nact_*nact_*nact_);
      SMITH::sort_indices<0,1,3,2,1,1,1,1>(four->data(), out->data(), nact_, nact_, nact_, nact_);
      // +/- (h)rdm3(i,j,x,k,y,z) * (<lx|yz>*)
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(rdm3->data(), six->data(), nact_*nact_, nact_, nact_, nact_*nact_); // (i,j,k,x,y,w)
      dgemm_("N", "T", nact_*nact_*nact_, nact_, nact_*nact_*nact_, sign, six->data(), nact_*nact_*nact_, ints->data(), nact_, 1.0, out->data(), nact_*nact_*nact_);
      return out;
    };

    kmat2_  = compute_kmat(nact_,  rdm2_,  rdm3_, fockact_c_, ints2_,  1.0);
    kmatp2_ = compute_kmat(nact_, hrdm2_, hrdm3_, fockact_h_, ints2_, -1.0);
  }
}


void NEVPT2::compute_abcd() {
  auto id2 = [this](                          const int k, const int l) { return         (        (k+nact_*l)); };
  auto id3 = [this](             const int j, const int k, const int l) { return         (j+nact_*(k+nact_*l)); };
  auto id4 = [this](const int i, const int j, const int k, const int l) { return i+nact_*(j+nact_*(k+nact_*l)); };
  // A matrices
  {
    shared_ptr<Matrix> amat2 = rdm2_->clone();
    {
      for (int b = 0; b != nact_; ++b)
        for (int a = 0; a != nact_; ++a)
          for (int bp = 0; bp != nact_; ++bp)
            for (int ap = 0; ap != nact_; ++ap)
              for (int c = 0; c != nact_; ++c) {
                amat2->element(ap+nact_*bp,a+nact_*b) += fockact_p_->element(c,a) * ardm2_->element(bp+nact_*ap,c+nact_*b) - fockact_p_->element(c,b) * ardm2_->element(bp+nact_*ap,a+nact_*c);
                for (int d = 0; d != nact_; ++d)
                  for (int e = 0; e != nact_; ++e)
                    amat2->element(ap+nact_*bp,a+nact_*b) += 0.5 * ints2_->element(c+nact_*d,e+nact_*a) * (ardm3_->element(id3(bp,ap,c),id3(e,d,b))
                                                                                                    + ardm3_->element(id3(bp,ap,d),id3(b,c,e)))
                                                         - 0.5 * ints2_->element(b+nact_*c,e+nact_*d) * (ardm3_->element(id3(bp,ap,a),id3(e,c,d))
                                                                                                    + ardm3_->element(id3(bp,ap,c),id3(d,a,e)));
              }
      shared_ptr<Matrix> tmp = amat2->copy();
      SMITH::sort_indices<1,0,3,2,0,1,1,1>(tmp->data(), amat2->data(), nact_, nact_, nact_, nact_);
    }
    shared_ptr<Matrix> amat3 = rdm3_->clone();
    shared_ptr<Matrix> amat3t = rdm3_->clone();
    {
      for (int c = 0; c != nact_; ++c)
        for (int b = 0; b != nact_; ++b)
          for (int a = 0; a != nact_; ++a)
            for (int cp = 0; cp != nact_; ++cp)
              for (int bp = 0; bp != nact_; ++bp)
                for (int ap = 0; ap != nact_; ++ap)
                  for (int d = 0; d != nact_; ++d) {
                    amat3->element(id3(ap,bp,cp),id3(a,b,c)) += fockact_p_->element(d,a)*ardm3_->element(id3(cp,ap,bp),id3(b,d,c))
                                                              - fockact_p_->element(c,d)*ardm3_->element(id3(cp,ap,bp),id3(b,a,d))
                                                              - fockact_p_->element(b,d)*ardm3_->element(id3(cp,ap,bp),id3(d,a,c));
                    amat3t->element(id3(ap,bp,cp),id3(a,b,c))+= fockact_p_->element(d,a)*srdm3_->element(id3(cp,ap,bp),id3(b,d,c))
                                                              - fockact_p_->element(c,d)*srdm3_->element(id3(cp,ap,bp),id3(b,a,d))
                                                              + fockact_p_->element(b,d)*srdm3_->element(id3(cp,ap,bp),id3(d,a,c));
                    for (int e = 0; e != nact_; ++e) {
                      amat3->element(id3(ap,bp,cp),id3(a,b,c)) += ints2_->element(id2(c,d),id2(e,a))*ardm3_->element(id3(cp,ap,bp),id3(b,d,e))
                                                            - 0.5*ints2_->element(id2(d,e),id2(e,a))*ardm3_->element(id3(cp,ap,bp),id3(b,d,c))
                                                            - 0.5*ints2_->element(id2(c,d),id2(d,e))*ardm3_->element(id3(cp,ap,bp),id3(b,a,e))
                                                            + 0.5*ints2_->element(id2(b,d),id2(d,e))*ardm3_->element(id3(cp,ap,bp),id3(e,a,c));
                      amat3t->element(id3(ap,bp,cp),id3(a,b,c))+= ints2_->element(id2(c,d),id2(e,a))*srdm3_->element(id3(cp,ap,bp),id3(b,d,e))
                                                            - 0.5*ints2_->element(id2(d,e),id2(e,a))*srdm3_->element(id3(cp,ap,bp),id3(b,d,c))
                                                            - 0.5*ints2_->element(id2(c,d),id2(d,e))*srdm3_->element(id3(cp,ap,bp),id3(b,a,e))
                                                            + 0.5*ints2_->element(id2(b,d),id2(d,e))*srdm3_->element(id3(cp,ap,bp),id3(e,a,c));
                      for (int f = 0; f != nact_; ++f) {
                        amat3->element(id3(ap,bp,cp),id3(a,b,c)) += ints2_->element(id2(d,e),id2(f,a))*ardm4_->element(id4(cp,ap,bp,b),id4(d,f,e,c))
                                                                  - ints2_->element(id2(d,c),id2(f,e))*ardm4_->element(id4(cp,ap,bp,b),id4(d,f,a,e))
                                                                  - ints2_->element(id2(d,b),id2(f,e))*ardm4_->element(id4(cp,ap,bp,e),id4(d,f,a,c));
                        amat3t->element(id3(ap,bp,cp),id3(a,b,c))+= ints2_->element(id2(d,e),id2(f,a))
                                                                          *((b == bp ? 2.0 : 0.0)*ardm3_->element(id3(cp,ap,d),id3(f,e,c)) - ardm4_->element(id4(cp,ap,b,bp),id4(d,f,e,c)))
                                                                  - ints2_->element(id2(d,c),id2(f,e))
                                                                          *((b == bp ? 2.0 : 0.0)*ardm3_->element(id3(cp,ap,d),id3(f,a,e)) - ardm4_->element(id4(cp,ap,b,bp),id4(d,f,a,e)))
                                                                  + ints2_->element(id2(d,b),id2(f,e))
                                                                          *((bp == e ? 2.0 : 0.0)*ardm3_->element(id3(cp,ap,d),id3(f,a,c)) - ardm4_->element(id4(cp,ap,e,bp),id4(d,f,a,c)));
                      }
                    }
                  }
    }
    amat2_ = amat2;
    amat3_ = amat3;
    amat3t_ = amat3t;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // B matrices
  {
    shared_ptr<Matrix> bmat2 = make_shared<Matrix>(nact_*nact_*nact_, nact_);
    shared_ptr<Matrix> bmat2t = make_shared<Matrix>(nact_*nact_*nact_, nact_);
    for (int a = 0; a != nact_; ++a)
      for (int cp = 0; cp != nact_; ++cp)
        for (int bp = 0; bp != nact_; ++bp)
          for (int ap = 0; ap != nact_; ++ap)
            for (int c = 0; c != nact_; ++c) {
              bmat2->element(id3(ap,bp,cp),a) -= (fockact_p_->element(a,c)*2.0-fockact_c_->element(a,c)) * ardm2_->element(id2(cp,ap),id2(bp,c));
              bmat2t->element(id3(ap,bp,cp),a) += fockact_c_->element(c,a) * srdm2_->element(id2(cp,ap),id2(bp,c));
              for (int e = 0; e != nact_; ++e)
                for (int f = 0; f != nact_; ++f) {
                  bmat2->element(id3(ap,bp,cp),a) -= ints2_->element(id2(a,c),id2(e,f)) * ardm3_->element(id3(cp,ap,bp),id3(e,c,f));
                  bmat2t->element(id3(ap,bp,cp),a) += ints2_->element(id2(a,c),id2(e,f)) * srdm3_->element(id3(cp,ap,bp),id3(e,c,f));
                }
            }
    bmat2_ = bmat2;
    bmat2t_ = bmat2t;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // C matrices
  {
    shared_ptr<Matrix> cmat2 = make_shared<Matrix>(nact_, nact_*nact_*nact_);
    shared_ptr<Matrix> cmat2t = make_shared<Matrix>(nact_, nact_*nact_*nact_);
    // <d c+ b+ a>
    shared_ptr<Matrix> s2rdm2 = ardm2_->clone();
    for (int a = 0; a != nact_; ++a)
      for (int b = 0; b != nact_; ++b)
        for (int c = 0; c != nact_; ++c)
          for (int d = 0; d != nact_; ++d)
            s2rdm2->element(id2(d,c),id2(b,a)) += (d == c ? 2.0 : 0.0) * rdm1_->element(b,a) - ardm2_->element(id2(c,d), id2(b,a));

    for (int c = 0; c != nact_; ++c)
      for (int a = 0; a != nact_; ++a)
        for (int b = 0; b != nact_; ++b)
          for (int ap = 0; ap != nact_; ++ap)
            for (int d = 0; d != nact_; ++d) {
              cmat2->element(ap, id3(a,b,c)) += fockact_p_->element(d,a)*ardm2_->element(id2(ap,b),id2(d,c)) - fockact_p_->element(c,d)*ardm2_->element(id2(ap,b),id2(a,d))
                                              - fockact_p_->element(b,d)*ardm2_->element(id2(ap,d),id2(a,c));
              cmat2t->element(ap, id3(a,b,c)) += fockact_p_->element(d,a)*s2rdm2->element(id2(ap,b),id2(d,c)) - fockact_p_->element(c,d)*s2rdm2->element(id2(ap,b),id2(a,d))
                                               + fockact_p_->element(b,d)*s2rdm2->element(id2(ap,d),id2(a,c));
              for (int e = 0; e != nact_; ++e) {
                for (int f = 0; f != nact_; ++f) {
                  cmat2->element(ap, id3(a,b,c)) += ints2_->element(id2(d,e),id2(f,a))*ardm3_->element(id3(ap,b,d),id3(f,e,c))
                                                  - ints2_->element(id2(d,c),id2(f,e))*ardm3_->element(id3(ap,b,d),id3(f,a,e))
                                                  - ints2_->element(id2(d,b),id2(f,e))*ardm3_->element(id3(ap,e,d),id3(f,a,c));
                  cmat2t->element(ap, id3(a,b,c))+= ints2_->element(id2(d,e),id2(f,a))*((ap == b ? 2.0 : 0.0)*ardm2_->element(id2(d,f),id2(e,c)) - ardm3_->element(id3(b,ap,d),id3(f,e,c)))
                                                  - ints2_->element(id2(d,c),id2(f,e))*((ap == b ? 2.0 : 0.0)*ardm2_->element(id2(d,f),id2(a,e)) - ardm3_->element(id3(b,ap,d),id3(f,a,e)))
                                                  + ints2_->element(id2(d,b),id2(f,e))*((ap == e ? 2.0 : 0.0)*ardm2_->element(id2(d,f),id2(a,c)) - ardm3_->element(id3(e,ap,d),id3(f,a,c)));
                }
                cmat2->element(ap, id3(a,b,c)) += ints2_->element(id2(c,d),id2(e,a))*ardm2_->element(id2(ap,b),id2(d,e))
                                             -0.5*ints2_->element(id2(d,e),id2(e,a))*ardm2_->element(id2(ap,b),id2(d,c))
                                             -0.5*ints2_->element(id2(c,d),id2(d,e))*ardm2_->element(id2(ap,b),id2(a,e))
                                             +0.5*ints2_->element(id2(b,d),id2(d,e))*ardm2_->element(id2(ap,e),id2(a,c));
                cmat2t->element(ap, id3(a,b,c))+= ints2_->element(id2(c,d),id2(e,a))*s2rdm2->element(id2(ap,b),id2(d,e))
                                             -0.5*ints2_->element(id2(d,e),id2(e,a))*s2rdm2->element(id2(ap,b),id2(d,c))
                                             -0.5*ints2_->element(id2(c,d),id2(d,e))*s2rdm2->element(id2(ap,b),id2(a,e))
                                             +0.5*ints2_->element(id2(b,d),id2(d,e))*s2rdm2->element(id2(ap,e),id2(a,c));
              }
            }
    cmat2_ = cmat2;
    cmat2t_ = cmat2t;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // D matrices
  {
    shared_ptr<Matrix> dmat2 = rdm2_->clone();
    for (int b = 0; b != nact_; ++b)
      for (int a = 0; a != nact_; ++a)
        for (int bp = 0; bp != nact_; ++bp)
          for (int ap = 0; ap != nact_; ++ap)
            for (int c = 0; c != nact_; ++c) {
              dmat2->element(ap+nact_*bp,a+nact_*b) -= fockact_p_->element(b,c) * ((a == ap ? 2.0 : 0.0) * rdm1_->element(bp,c) - rdm2_->element(a+nact_*bp,ap+nact_*c));
              dmat2->element(ap+nact_*bp,a+nact_*b) += fockact_p_->element(a,c) * ((c == ap ? 2.0 : 0.0) * rdm1_->element(bp,b) - rdm2_->element(c+nact_*bp,ap+nact_*b));
              for (int d = 0; d != nact_; ++d)
                for (int e = 0; e != nact_; ++e) {
                  dmat2->element(ap+nact_*bp,a+nact_*b) -= 0.5 * ints2_->element(c+nact_*b,e+nact_*d)
                           * ((a == ap ? 2.0 : 0.0) * ardm2_->element(c+nact_*e,bp+nact_*d) - ardm3_->element(id3(c,e,a),id3(ap,bp,d))
                           +  (a == ap ? 2.0 : 0.0) * ardm2_->element(bp+nact_*d,c+nact_*e) - ardm3_->element(id3(a,ap,bp),id3(d,c,e))
                           + (ap == bp ? 1.0 : 0.0) *(ardm2_->element(c+nact_*e,a+nact_*d)  + ardm2_->element(a+nact_*d,c+nact_*e))
                           +  (c == ap ? 1.0 : 0.0) *((a == e  ? 2.0 : 0.0) * rdm1_->element(bp, d) - ardm2_->element(a+nact_*e,bp+nact_*d))
                           - (bp == e  ? 1.0 : 0.0) *((a == ap ? 2.0 : 0.0) * rdm1_->element(c,d)   - ardm2_->element(a+nact_*ap,c+nact_*d)));
                  dmat2->element(ap+nact_*bp,a+nact_*b) += 0.5 * ints2_->element(c+nact_*d,e+nact_*a)
                           * ((d == ap ? 2.0 : 0.0) * ardm2_->element(c+nact_*e,bp+nact_*b) - ardm3_->element(id3(c,e,d),id3(ap,bp,b))
                           +  (d == ap ? 2.0 : 0.0) * ardm2_->element(bp+nact_*b,c+nact_*e) - ardm3_->element(id3(d,ap,bp),id3(b,c,e))
                           + (ap == bp ? 1.0 : 0.0) *(ardm2_->element(c+nact_*e,d+nact_*b)  + ardm2_->element(d+nact_*b,c+nact_*e))
                           +  (c == ap ? 1.0 : 0.0) *((d == e  ? 2.0 : 0.0) * rdm1_->element(bp, b) - ardm2_->element(d+nact_*e,bp+nact_*b))
                           - (bp == e  ? 1.0 : 0.0) *((d == ap ? 2.0 : 0.0) * rdm1_->element(c,b)   - ardm2_->element(d+nact_*ap,c+nact_*b)));
                }
            }
    shared_ptr<Matrix> tmp = dmat2->copy();
    SMITH::sort_indices<1,0,3,2,0,1,1,1>(tmp->data(), dmat2->data(), nact_, nact_, nact_, nact_);
    dmat2_ = dmat2;
  }
  {
    shared_ptr<Matrix> dmat1 = rdm1_->clone();
    shared_ptr<Matrix> dmat1t = rdm1_->clone();
    for (int a = 0; a != nact_; ++a)
      for (int ap = 0; ap != nact_; ++ap)
        for (int c = 0; c != nact_; ++c) {
          dmat1->element(ap,a) += -(fockact_p_->element(a,c)*2.0-fockact_c_->element(a,c)) * rdm1_->element(c,ap);
          dmat1t->element(ap,a) += fockact_c_->element(a,c) * hrdm1_->element(c,ap);
          for (int e = 0; e != nact_; ++e)
            for (int f = 0; f != nact_; ++f) {
              dmat1->element(ap,a) += - ints2_->element(id2(a,c),id2(e,f)) * ardm2_->element(id2(ap,e),id2(c,f));
              dmat1t->element(ap,a) += ints2_->element(id2(a,c),id2(e,f)) * ((ap == e ? 2.0 : 0.0)*rdm1_->element(c,f) - ardm2_->element(id2(e,ap),id2(c,f)));
            }
        }
    dmat1_ = dmat1;
    dmat1t_ = dmat1t;
  }
}

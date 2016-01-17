//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: nevpt2_mat.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifdef NEVPT2IMPL

#ifndef NDEBUG
template<typename T> inline bool is_hermitian_conjugate(shared_ptr<const T> t, shared_ptr<const T> u) { assert(false); return false; }
template<> inline bool is_hermitian_conjugate<Matrix>(shared_ptr<const Matrix> t, shared_ptr<const Matrix> u) {
  return abs((*t->transpose()-*u).rms()) < 1.0e-8;
}
template<> inline bool is_hermitian_conjugate<ZMatrix>(shared_ptr<const ZMatrix> t, shared_ptr<const ZMatrix> u) {
  return abs((*t->transpose_conjg()-*u).rms()) < 1.0e-8;
}
#endif

template<typename DataType>
void NEVPT2<DataType>::compute_kmat() {
  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;
  {
    // Eq. (27)
    auto kmat = make_shared<MatType>(nact_, nact_, true);
    btas::contract(1.0, *fockact_c_, {0,1}, *rdm1_, {2,1}, 0.0, *kmat, {0,2});
    *kmat += *qvec_;
    // transpose (assuming that kmat is hermitian...)
    kmat = kmat->get_conjg();
    kmat->localize();

    // Eq. (A4)
    auto kmatp = make_shared<MatType>(*kmat * (-1.0));
    *kmatp += *fockact_->get_conjg() * fac2;
    kmatp->localize();

    kmat_ = kmat;
    kmatp_ = kmatp;
  }
  {
    auto compute_kmat = [this,&fac2](shared_ptr<const MatType> rdm2, shared_ptr<const MatType> rdm3, shared_ptr<const MatType> fock, const double sign) {
      auto out = rdm2->clone();
      // Eq. (A7) and (A9)
      for (int b = 0; b != nact_; ++b)
        for (int a = 0; a != nact_; ++a)
          for (int bp = 0; bp != nact_; ++bp)
            for (int ap = 0; ap != nact_; ++ap)
              for (int c = 0; c != nact_; ++c) {
                out->element(ap+nact_*bp, a+nact_*b) += rdm2->element(ap+nact_*bp, c+nact_*b) * fock->element(a, c)
                                                      + rdm2->element(ap+nact_*bp, a+nact_*c) * fock->element(b, c);
                for (int d = 0; d != nact_; ++d)
                  for (int e = 0; e != nact_; ++e) {
                    out->element(ap+nact_*bp, a+nact_*b) += 0.5 * ints2_->element(e+nact_*a, c+nact_*d)
                                                           * (sign*2.0 * rdm3->element(ap+nact_*(bp+nact_*e), d+nact_*(b+nact_*c))
                                                             +sign*(b == e ? rdm2->element(ap+nact_*bp, d+nact_*c) : 0.0))
                                                          + 0.5 * ints2_->element(e+nact_*b, c+nact_*d)
                                                           * (sign*2.0 * rdm3->element(ap+nact_*(bp+nact_*e), a+nact_*(d+nact_*c))
                                                             +sign*(a == e ? rdm2->element(ap+nact_*bp, c+nact_*d) : 0.0));
                  }
              }
      return out;
    };

    kmat2_  = compute_kmat( rdm2_,  rdm3_, fockact_c_,  1.0);
    kmatp2_ = compute_kmat(hrdm2_, hrdm3_, fockact_h_, -1.0);
    // for CAS references, these are Hermitian; if not, suspect the factor of 2.0 above
    assert(kmat_->is_hermitian());
    assert(kmatp_->is_hermitian());
    assert(kmat2_->is_hermitian());
    assert(kmatp2_->is_hermitian());
  }
}


template<typename DataType>
void NEVPT2<DataType>::compute_abcd() {
  auto id2 = [this](                          const int k, const int l) { return         (        (k+nact_*l)); };
  auto id3 = [this](             const int j, const int k, const int l) { return         (j+nact_*(k+nact_*l)); };
  auto id4 = [this](const int i, const int j, const int k, const int l) { return i+nact_*(j+nact_*(k+nact_*l)); };

  const double fac2 = is_same<DataType,double>::value ? 2.0 : 1.0;
  // A matrices
  {
    shared_ptr<MatType> amat2 = rdm2_->clone();
    {
      for (int b = 0; b != nact_; ++b)
        for (int a = 0; a != nact_; ++a)
          for (int bp = 0; bp != nact_; ++bp)
            for (int ap = 0; ap != nact_; ++ap)
              for (int c = 0; c != nact_; ++c) {
                amat2->element(ap+nact_*bp,a+nact_*b) += fockact_p_->element(c,a) * ardm2_->element(bp+nact_*ap,c+nact_*b)
                                                       - fockact_p_->element(b,c) * ardm2_->element(bp+nact_*ap,a+nact_*c);
                for (int d = 0; d != nact_; ++d)
                  for (int e = 0; e != nact_; ++e)
                    amat2->element(ap+nact_*bp,a+nact_*b) += 0.5 * ints2_->element(c+nact_*d,e+nact_*a) * (ardm3_->element(id3(bp,ap,c),id3(e,d,b))
                                                                                                         + ardm3_->element(id3(bp,ap,d),id3(b,c,e)))
                                                           - 0.5 * ints2_->element(b+nact_*c,e+nact_*d) * (ardm3_->element(id3(bp,ap,a),id3(e,c,d))
                                                                                                         + ardm3_->element(id3(bp,ap,c),id3(d,a,e)));
              }
      shared_ptr<MatType> tmp = amat2->copy();
      sort_indices<1,0,3,2,0,1,1,1>(tmp->data(), amat2->data(), nact_, nact_, nact_, nact_);
    }
    shared_ptr<MatType> amat3 = rdm3_->clone();
    shared_ptr<MatType> amat3t = rdm3_->clone();
    {
      const MatType fock_pp = *fockact_p_ * 2.0 - *fockact_c_;
      for (int c = 0; c != nact_; ++c)
        for (int b = 0; b != nact_; ++b)
          for (int a = 0; a != nact_; ++a)
            for (int cp = 0; cp != nact_; ++cp)
              for (int bp = 0; bp != nact_; ++bp)
                for (int ap = 0; ap != nact_; ++ap)
                  for (int d = 0; d != nact_; ++d) {
                    amat3->element(id3(ap,bp,cp),id3(a,b,c)) += fock_pp.element(d,a)*ardm3_->element(id3(cp,ap,bp),id3(b,d,c))
                                                              - fockact_c_->element(c,d)*ardm3_->element(id3(cp,ap,bp),id3(b,a,d))
                                                              - fock_pp.element(b,d)*ardm3_->element(id3(cp,ap,bp),id3(d,a,c));
                    amat3t->element(id3(ap,bp,cp),id3(a,b,c))+= fock_pp.element(d,a)*srdm3_->element(id3(cp,ap,bp),id3(b,d,c))
                                                              - fockact_c_->element(c,d)*srdm3_->element(id3(cp,ap,bp),id3(b,a,d))
                                                              + fockact_c_->element(d,b)*srdm3_->element(id3(cp,ap,bp),id3(d,a,c));
                    for (int e = 0; e != nact_; ++e) {
                      amat3->element(id3(ap,bp,cp),id3(a,b,c)) += ints2_->element(id2(c,d),id2(e,a))*ardm3_->element(id3(cp,ap,bp),id3(b,d,e));
                      amat3t->element(id3(ap,bp,cp),id3(a,b,c))+= ints2_->element(id2(c,d),id2(e,a))*srdm3_->element(id3(cp,ap,bp),id3(b,d,e));
                      for (int f = 0; f != nact_; ++f) {
                        amat3->element(id3(ap,bp,cp),id3(a,b,c)) += ints2_->element(id2(d,e),id2(f,a))*ardm4_->element(id4(cp,ap,bp,b),id4(d,f,e,c))
                                                                  - ints2_->element(id2(d,c),id2(f,e))*ardm4_->element(id4(cp,ap,bp,b),id4(d,f,a,e))
                                                                  - ints2_->element(id2(d,b),id2(f,e))*ardm4_->element(id4(cp,ap,bp,e),id4(d,f,a,c));
                        amat3t->element(id3(ap,bp,cp),id3(a,b,c))+= ints2_->element(id2(d,e),id2(f,a))
                                                                          *((b == bp ? fac2 : 0.0)*ardm3_->element(id3(cp,ap,d),id3(f,e,c)) - ardm4_->element(id4(cp,ap,b,bp),id4(d,f,e,c)))
                                                                  - ints2_->element(id2(d,c),id2(f,e))
                                                                          *((b == bp ? fac2 : 0.0)*ardm3_->element(id3(cp,ap,d),id3(f,a,e)) - ardm4_->element(id4(cp,ap,b,bp),id4(d,f,a,e)))
                                                                  + ints2_->element(id2(d,e),id2(f,b))
                                                                          *((bp == e ? fac2 : 0.0)*ardm3_->element(id3(cp,ap,d),id3(f,a,c)) - ardm4_->element(id4(cp,ap,e,bp),id4(d,f,a,c)));
                      }
                    }
                  }
    }
    amat2_ = amat2;
    assert(amat2_->is_hermitian());
    amat3_ = amat3;
    assert(amat3_->is_hermitian());
    amat3t_ = amat3t;
    assert(amat3t_->is_hermitian());
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // B matrices
  {
    auto bmat2 = make_shared<MatType>(nact_*nact_*nact_, nact_, true);
    auto bmat2t = make_shared<MatType>(nact_*nact_*nact_, nact_, true);
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
                  bmat2t->element(id3(ap,bp,cp),a) += ints2_->element(id2(e,c),id2(a,f)) * srdm3_->element(id3(cp,ap,bp),id3(e,c,f));
                }
            }
    bmat2_ = bmat2;
    bmat2t_ = bmat2t;
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // C matrices
  {
    auto cmat2 = make_shared<MatType>(nact_, nact_*nact_*nact_, true);
    auto cmat2t = make_shared<MatType>(nact_, nact_*nact_*nact_, true);
    // <d c+ b+ a>
    shared_ptr<MatType> s2rdm2 = ardm2_->clone();
    for (int a = 0; a != nact_; ++a)
      for (int b = 0; b != nact_; ++b)
        for (int c = 0; c != nact_; ++c)
          for (int d = 0; d != nact_; ++d)
            s2rdm2->element(id2(d,c),id2(b,a)) += (d == c ? fac2 : 0.0) * rdm1_->element(b,a) - ardm2_->element(id2(c,d), id2(b,a));

    const MatType fock_pp = *fockact_p_ * 2.0 - *fockact_c_;
    for (int c = 0; c != nact_; ++c)
      for (int a = 0; a != nact_; ++a)
        for (int b = 0; b != nact_; ++b)
          for (int ap = 0; ap != nact_; ++ap)
            for (int d = 0; d != nact_; ++d) {
              cmat2->element(ap, id3(a,b,c)) += fock_pp.element(d,a)*ardm2_->element(id2(ap,b),id2(d,c)) - fockact_c_->element(c,d)*ardm2_->element(id2(ap,b),id2(a,d))
                                              - fock_pp.element(b,d)*ardm2_->element(id2(ap,d),id2(a,c));
              cmat2t->element(ap, id3(a,b,c)) += fock_pp.element(d,a)*s2rdm2->element(id2(ap,b),id2(d,c)) - fockact_c_->element(c,d)*s2rdm2->element(id2(ap,b),id2(a,d))
                                               + fockact_c_->element(d,b)*s2rdm2->element(id2(ap,d),id2(a,c));
              for (int e = 0; e != nact_; ++e) {
                for (int f = 0; f != nact_; ++f) {
                  cmat2->element(ap, id3(a,b,c)) += ints2_->element(id2(d,e),id2(f,a))*ardm3_->element(id3(ap,b,d),id3(f,e,c))
                                                  - ints2_->element(id2(d,c),id2(f,e))*ardm3_->element(id3(ap,b,d),id3(f,a,e))
                                                  - ints2_->element(id2(d,b),id2(f,e))*ardm3_->element(id3(ap,e,d),id3(f,a,c));
                  cmat2t->element(ap, id3(a,b,c))+= ints2_->element(id2(d,e),id2(f,a))*((ap == b ? fac2 : 0.0)*ardm2_->element(id2(d,f),id2(e,c)) - ardm3_->element(id3(b,ap,d),id3(f,e,c)))
                                                  - ints2_->element(id2(d,c),id2(f,e))*((ap == b ? fac2 : 0.0)*ardm2_->element(id2(d,f),id2(a,e)) - ardm3_->element(id3(b,ap,d),id3(f,a,e)))
                                                  + ints2_->element(id2(d,e),id2(f,b))*((ap == e ? fac2 : 0.0)*ardm2_->element(id2(d,f),id2(a,c)) - ardm3_->element(id3(e,ap,d),id3(f,a,c)));
                }
                cmat2->element(ap, id3(a,b,c)) += ints2_->element(id2(c,d),id2(e,a))*ardm2_->element(id2(ap,b),id2(d,e));
                cmat2t->element(ap, id3(a,b,c))+= ints2_->element(id2(c,d),id2(e,a))*s2rdm2->element(id2(ap,b),id2(d,e));
              }
            }
    cmat2_ = cmat2;
    cmat2t_ = cmat2t;
  }

  // checking
  assert(is_hermitian_conjugate(bmat2_, cmat2_));
  assert(is_hermitian_conjugate(bmat2t_, cmat2t_));

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // D matrices
  {
    shared_ptr<MatType> s2rdm2 = ardm2_->clone();
    shared_ptr<MatType> s2rdm3 = ardm3_->clone();
    for (int a = 0; a != nact_; ++a)
      for (int b = 0; b != nact_; ++b)
        for (int c = 0; c != nact_; ++c)
          for (int d = 0; d != nact_; ++d)
            s2rdm2->element(id2(d,c),id2(b,a)) += (d == c ? fac2 : 0.0) * rdm1_->element(b,a) - ardm2_->element(id2(c,d), id2(b,a));
    for (int a = 0; a != nact_; ++a)
      for (int b = 0; b != nact_; ++b)
        for (int c = 0; c != nact_; ++c)
          for (int d = 0; d != nact_; ++d)
            for (int e = 0; e != nact_; ++e)
              for (int f = 0; f != nact_; ++f)
                s2rdm3->element(id3(f,e,d),id3(c,b,a)) += (f == e ? fac2 : 0.0) * ardm2_->element(id2(d,c),id2(b,a)) - ardm3_->element(id3(e,f,d), id3(c,b,a));

    shared_ptr<MatType> dmat2 = rdm2_->clone();
    for (int b = 0; b != nact_; ++b)
      for (int a = 0; a != nact_; ++a)
        for (int bp = 0; bp != nact_; ++bp)
          for (int ap = 0; ap != nact_; ++ap)
            for (int c = 0; c != nact_; ++c) {
              dmat2->element(ap+nact_*bp,a+nact_*b) -= fockact_p_->element(b,c) * ((a == ap ? fac2: 0.0) * rdm1_->element(bp,c) - rdm2_->element(a+nact_*bp,ap+nact_*c));
              dmat2->element(ap+nact_*bp,a+nact_*b) += fockact_p_->element(c,a) * ((c == ap ? fac2: 0.0) * rdm1_->element(bp,b) - rdm2_->element(c+nact_*bp,ap+nact_*b));
              for (int d = 0; d != nact_; ++d)
                for (int e = 0; e != nact_; ++e) {
                  dmat2->element(ap+nact_*bp,a+nact_*b) -= 0.5 * ints2_->element(c+nact_*b,e+nact_*d)
                           * (srdm3_->element(id3(c,e,ap),id3(a,bp,d)) +  s2rdm3->element(id3(ap,a,bp),id3(d,c,e))
                           + (ap == bp ? 1.0 : 0.0) *(ardm2_->element(c+nact_*e,a+nact_*d)  + ardm2_->element(a+nact_*d,c+nact_*e))
                           +  (c == ap ? 1.0 : 0.0) *s2rdm2->element(e+nact_*a,bp+nact_*d)
                           - (bp == e  ? 1.0 : 0.0) *s2rdm2->element(ap+nact_*a,c+nact_*d));
                  dmat2->element(ap+nact_*bp,a+nact_*b) += 0.5 * ints2_->element(c+nact_*d,e+nact_*a)
                           * (srdm3_->element(id3(c,e,ap),id3(d,bp,b)) +  s2rdm3->element(id3(ap,d,bp),id3(b,c,e))
                           + (ap == bp ? 1.0 : 0.0) *(ardm2_->element(c+nact_*e,d+nact_*b)  + ardm2_->element(d+nact_*b,c+nact_*e))
                           +  (c == ap ? 1.0 : 0.0) * s2rdm2->element(id2(e,d), id2(bp,b))
                           - (bp == e  ? 1.0 : 0.0) * s2rdm2->element(id2(ap,d), id2(c,b)));

                }
            }
    shared_ptr<MatType> tmp = dmat2->copy();
    sort_indices<1,0,3,2,0,1,1,1>(tmp->data(), dmat2->data(), nact_, nact_, nact_, nact_);
    dmat2_ = dmat2;
    assert(dmat2_->is_hermitian());
  }
  {
    shared_ptr<MatType> dmat1 = rdm1_->clone();
    shared_ptr<MatType> dmat1t = rdm1_->clone();
    for (int a = 0; a != nact_; ++a)
      for (int ap = 0; ap != nact_; ++ap)
        for (int c = 0; c != nact_; ++c) {
          dmat1->element(ap,a) += -(fockact_p_->element(a,c)*2.0-fockact_c_->element(a,c)) * rdm1_->element(ap,c);
          dmat1t->element(ap,a) += fockact_c_->element(c,a) * hrdm1_->element(c,ap);
          for (int e = 0; e != nact_; ++e)
            for (int f = 0; f != nact_; ++f) {
              dmat1->element(ap,a) += - ints2_->element(id2(a,c),id2(e,f)) * ardm2_->element(id2(ap,e),id2(c,f));
              dmat1t->element(ap,a) += ints2_->element(id2(e,c),id2(a,f)) * ((ap == e ? fac2 : 0.0)*rdm1_->element(c,f) - ardm2_->element(id2(e,ap),id2(c,f)));
            }
        }
    dmat1_ = dmat1;
    assert(dmat1_->is_hermitian());
    dmat1t_ = dmat1t;
    assert(dmat1t_->is_hermitian());
  }
}

#endif

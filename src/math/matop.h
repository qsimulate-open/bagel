//
// BAGEL - Parallel electron correlation program.
// Filename: matop.h
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

#ifndef __SRC_MATH_MATOP_H
#define __SRC_MATH_MATOP_H

#include <src/math/vectorb.h>
#include <src/math/matrix.h>
#include <src/math/zmatrix.h>

namespace bagel {

// operator+= and -=
inline Matrix&  operator+=(Matrix& a,  const Matrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline Matrix&  operator+=(Matrix& a,  const MatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline MatView& operator+=(MatView& a, const Matrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline MatView& operator+=(MatView& a, const MatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline Matrix&  operator-=(Matrix& a,  const Matrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline Matrix&  operator-=(Matrix& a,  const MatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline MatView& operator-=(MatView& a, const Matrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline MatView& operator-=(MatView& a, const MatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline VectorB& operator+=(VectorB& a, const VectorB& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline VectorB& operator+=(VectorB& a, const VecView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline VecView& operator+=(VecView& a, const VectorB& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline VecView& operator+=(VecView& a, const VecView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline VectorB& operator-=(VectorB& a, const VectorB& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline VectorB& operator-=(VectorB& a, const VecView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline VecView& operator-=(VecView& a, const VectorB& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline VecView& operator-=(VecView& a, const VecView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }

inline ZMatrix&  operator+=(ZMatrix& a,  const ZMatrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatrix&  operator+=(ZMatrix& a,  const ZMatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatView& operator+=(ZMatView& a, const ZMatrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatView& operator+=(ZMatView& a, const ZMatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatrix&  operator-=(ZMatrix& a,  const ZMatrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatrix&  operator-=(ZMatrix& a,  const ZMatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatView& operator-=(ZMatView& a, const ZMatrix& b)  { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZMatView& operator-=(ZMatView& a, const ZMatView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZVectorB& operator+=(ZVectorB& a, const ZVectorB& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZVectorB& operator+=(ZVectorB& a, const ZVecView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZVecView& operator+=(ZVecView& a, const ZVectorB& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZVecView& operator+=(ZVecView& a, const ZVecView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }
inline ZVectorB& operator-=(ZVectorB& a, const ZVectorB& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZVectorB& operator-=(ZVectorB& a, const ZVecView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZVecView& operator-=(ZVecView& a, const ZVectorB& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }
inline ZVecView& operator-=(ZVecView& a, const ZVecView& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }

// operator+ and -
inline Matrix operator+(const Matrix& a,  const Matrix& b)  { Matrix out(a); out += b; return out; }
inline Matrix operator+(const Matrix& a,  const MatView& b) { Matrix out(a); out += b; return out; }
inline Matrix operator+(const MatView& a, const Matrix& b)  { Matrix out(a); out += b; return out; }
inline Matrix operator+(const MatView& a, const MatView& b) { Matrix out(a); out += b; return out; }
inline Matrix operator-(const Matrix& a,  const Matrix& b)  { Matrix out(a); out -= b; return out; }
inline Matrix operator-(const Matrix& a,  const MatView& b) { Matrix out(a); out -= b; return out; }
inline Matrix operator-(const MatView& a, const Matrix& b)  { Matrix out(a); out -= b; return out; }
inline Matrix operator-(const MatView& a, const MatView& b) { Matrix out(a); out -= b; return out; }
inline VectorB operator+(const VectorB& a, const VectorB& b) { VectorB out(a); out += b; return out; }
inline VectorB operator+(const VectorB& a, const VecView& b) { VectorB out(a); out += b; return out; }
inline VectorB operator+(const VecView& a, const VectorB& b) { VectorB out(a); out += b; return out; }
inline VectorB operator+(const VecView& a, const VecView& b) { VectorB out(a); out += b; return out; }
inline VectorB operator-(const VectorB& a, const VectorB& b) { VectorB out(a); out -= b; return out; }
inline VectorB operator-(const VectorB& a, const VecView& b) { VectorB out(a); out -= b; return out; }
inline VectorB operator-(const VecView& a, const VectorB& b) { VectorB out(a); out -= b; return out; }
inline VectorB operator-(const VecView& a, const VecView& b) { VectorB out(a); out -= b; return out; }

inline ZMatrix operator+(const ZMatrix& a,  const ZMatrix& b)  { ZMatrix out(a); out += b; return out; }
inline ZMatrix operator+(const ZMatrix& a,  const ZMatView& b) { ZMatrix out(a); out += b; return out; }
inline ZMatrix operator+(const ZMatView& a, const ZMatrix& b)  { ZMatrix out(a); out += b; return out; }
inline ZMatrix operator+(const ZMatView& a, const ZMatView& b) { ZMatrix out(a); out += b; return out; }
inline ZMatrix operator-(const ZMatrix& a,  const ZMatrix& b)  { ZMatrix out(a); out -= b; return out; }
inline ZMatrix operator-(const ZMatrix& a,  const ZMatView& b) { ZMatrix out(a); out -= b; return out; }
inline ZMatrix operator-(const ZMatView& a, const ZMatrix& b)  { ZMatrix out(a); out -= b; return out; }
inline ZMatrix operator-(const ZMatView& a, const ZMatView& b) { ZMatrix out(a); out -= b; return out; }
inline ZVectorB operator+(const ZVectorB& a, const ZVectorB& b) { ZVectorB out(a); out += b; return out; }
inline ZVectorB operator+(const ZVectorB& a, const ZVecView& b) { ZVectorB out(a); out += b; return out; }
inline ZVectorB operator+(const ZVecView& a, const ZVectorB& b) { ZVectorB out(a); out += b; return out; }
inline ZVectorB operator+(const ZVecView& a, const ZVecView& b) { ZVectorB out(a); out += b; return out; }
inline ZVectorB operator-(const ZVectorB& a, const ZVectorB& b) { ZVectorB out(a); out -= b; return out; }
inline ZVectorB operator-(const ZVectorB& a, const ZVecView& b) { ZVectorB out(a); out -= b; return out; }
inline ZVectorB operator-(const ZVecView& a, const ZVectorB& b) { ZVectorB out(a); out -= b; return out; }
inline ZVectorB operator-(const ZVecView& a, const ZVecView& b) { ZVectorB out(a); out -= b; return out; }

namespace impl {

  template<class A, class B,
           class = typename std::enable_if<std::is_same<typename A::value_type, typename B::value_type>::value
                                       and std::is_same<typename A::value_type, double>::value
                                       and btas::is_boxtensor<A>::value and btas::is_boxtensor<B>::value
                                          >::type
          >
  Matrix multNN(const A& a, const B& b) {
    assert(a.rank() == 2 && b.rank() == 2);
    const int l = a.ndim();
    assert(a.mdim() == b.ndim());
    const int n = b.mdim();
    Matrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    const int m = a.mdim();
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      btas::contract(1.0, a, {0,2}, b, {2,1}, 0.0, out, {0,1});
#ifdef HAVE_SCALAPACK
    } else {
      std::unique_ptr<double[]> locala = a.getlocal();
      std::unique_ptr<double[]> localb = b.getlocal();
      std::unique_ptr<double[]> localc = out.getlocal();
      pdgemm_("N", "N", l, n, m, 1.0, locala.get(), a.desc().data(), localb.get(), b.desc().data(), 0.0, localc.get(), out.desc().data());
      out.setlocal(localc);
    }
#endif
    return out;
  }
  template<class A, class B,
           class = typename std::enable_if<std::is_same<typename A::value_type, typename B::value_type>::value
                                       and std::is_same<typename A::value_type, double>::value
                                       and btas::is_boxtensor<A>::value and btas::is_boxtensor<B>::value
                                          >::type
          >
  Matrix multTN(const A& a, const B& b) {
    assert(a.rank() == 2 && b.rank() == 2);
    const int l = a.mdim();
    assert(a.ndim() == b.ndim());
    const int n = b.mdim();
    Matrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    const int m = a.ndim();
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      btas::contract(1.0, a, {2,0}, b, {2,1}, 0.0, out, {0,1});
#ifdef HAVE_SCALAPACK
    } else {
      std::unique_ptr<double[]> locala = a.getlocal();
      std::unique_ptr<double[]> localb = b.getlocal();
      std::unique_ptr<double[]> localc = out.getlocal();
      pdgemm_("T", "N", l, n, m, 1.0, locala.get(), a.desc().data(), localb.get(), b.desc().data(), 0.0, localc.get(), out.desc().data());
      out.setlocal(localc);
    }
#endif

    return out;
  }

  template<class A, class B,
           class = typename std::enable_if<std::is_same<typename A::value_type, typename B::value_type>::value
                                       and std::is_same<typename A::value_type, double>::value
                                       and btas::is_boxtensor<A>::value and btas::is_boxtensor<B>::value
                                          >::type
          >
  Matrix multNT(const A& a, const B& b) {
    assert(a.rank() == 2 && b.rank() == 2);
    const int l = a.ndim();
    assert(a.mdim() == b.mdim());
    const int n = b.ndim();

    Matrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    const int m = a.mdim();
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      btas::contract(1.0, a, {0,2}, b, {1,2}, 0.0, out, {0,1});
#ifdef HAVE_SCALAPACK
    } else {
      std::unique_ptr<double[]> locala = a.getlocal();
      std::unique_ptr<double[]> localb = b.getlocal();
      std::unique_ptr<double[]> localc = out.getlocal();
      pdgemm_("N", "T", l, n, m, 1.0, locala.get(), a.desc().data(), localb.get(), b.desc().data(), 0.0, localc.get(), out.desc().data());
      out.setlocal(localc);
    }
#endif

    return out;
  }

  template<class A, class B,
           class = typename std::enable_if<std::is_same<typename A::value_type, typename B::value_type>::value
                                       and std::is_same<typename A::value_type, std::complex<double>>::value
                                       and btas::is_boxtensor<A>::value and btas::is_boxtensor<B>::value
                                          >::type
          >
  ZMatrix multNN(const A& a, const B& b) {
    assert(a.rank() == 2 && b.rank() == 2);
    const int l = a.ndim();
    const int m = a.mdim();
    assert(a.mdim() == b.ndim());
    const int n = b.mdim();
    ZMatrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      zgemm3m_("N", "N", l, n, m, 1.0, a.data(), l, b.data(), m, 0.0, out.data(), l);
#ifdef HAVE_SCALAPACK
    } else {
      std::unique_ptr<std::complex<double>[]> locala = a.getlocal();
      std::unique_ptr<std::complex<double>[]> localb = b.getlocal();
      std::unique_ptr<std::complex<double>[]> localc = out.getlocal();
      pzgemm_("N", "N", l, n, m, 1.0, locala.get(), a.desc().data(), localb.get(), b.desc().data(), 0.0, localc.get(), out.desc().data());
      out.setlocal(localc);
    }
#endif
    return out;
  }
  template<class A, class B,
           class = typename std::enable_if<std::is_same<typename A::value_type, typename B::value_type>::value
                                       and std::is_same<typename A::value_type, std::complex<double>>::value
                                       and btas::is_boxtensor<A>::value and btas::is_boxtensor<B>::value
                                          >::type
          >
  ZMatrix multTN(const A& a, const B& b) {
    assert(a.rank() == 2 && b.rank() == 2);
    const int l = a.mdim();
    const int m = a.ndim();
    assert(a.ndim() == b.ndim());
    const int n = b.mdim();
    ZMatrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      zgemm3m_("C", "N", l, n, m, 1.0, a.data(), m, b.data(), m, 0.0, out.data(), l);
#ifdef HAVE_SCALAPACK
    } else {
      std::unique_ptr<std::complex<double>[]> locala = a.getlocal();
      std::unique_ptr<std::complex<double>[]> localb = b.getlocal();
      std::unique_ptr<std::complex<double>[]> localc = out.getlocal();
      pzgemm_("C", "N", l, n, m, 1.0, locala.get(), a.desc().data(), localb.get(), b.desc().data(), 0.0, localc.get(), out.desc().data());
      out.setlocal(localc);
    }
#endif

    return out;
  }
  template<class A, class B,
           class = typename std::enable_if<std::is_same<typename A::value_type, typename B::value_type>::value
                                       and std::is_same<typename A::value_type, std::complex<double>>::value
                                       and btas::is_boxtensor<A>::value and btas::is_boxtensor<B>::value
                                          >::type
          >
  ZMatrix multNT(const A& a, const B& b) {
    assert(a.rank() == 2 && b.rank() == 2);
    const int l = a.ndim();
    const int m = a.mdim();
    assert(a.mdim() == b.mdim());
    const int n = b.ndim();
    ZMatrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      zgemm3m_("N", "C", l, n, m, 1.0, a.data(), l, b.data(), n, 0.0, out.data(), l);
#ifdef HAVE_SCALAPACK
    } else {
      std::unique_ptr<std::complex<double>[]> locala = a.getlocal();
      std::unique_ptr<std::complex<double>[]> localb = b.getlocal();
      std::unique_ptr<std::complex<double>[]> localc = out.getlocal();
      pzgemm_("N", "C", l, n, m, 1.0, locala.get(), a.desc().data(), localb.get(), b.desc().data(), 0.0, localc.get(), out.desc().data());
      out.setlocal(localc);
    }
#endif

    return out;
  }

}

// WARNING - we are abusing the operator overload!

// operator*
inline Matrix operator*(const Matrix& a,  const Matrix& b)  { return impl::multNN(a,b); }
inline Matrix operator*(const Matrix& a,  const MatView& b) { return impl::multNN(a,b); }
inline Matrix operator*(const MatView& a, const Matrix& b)  { return impl::multNN(a,b); }
inline Matrix operator*(const MatView& a, const MatView& b) { return impl::multNN(a,b); }
inline Matrix operator%(const Matrix& a,  const Matrix& b)  { return impl::multTN(a,b); }
inline Matrix operator%(const Matrix& a,  const MatView& b) { return impl::multTN(a,b); }
inline Matrix operator%(const MatView& a, const Matrix& b)  { return impl::multTN(a,b); }
inline Matrix operator%(const MatView& a, const MatView& b) { return impl::multTN(a,b); }
inline Matrix operator^(const Matrix& a,  const Matrix& b)  { return impl::multNT(a,b); }
inline Matrix operator^(const Matrix& a,  const MatView& b) { return impl::multNT(a,b); }
inline Matrix operator^(const MatView& a, const Matrix& b)  { return impl::multNT(a,b); }
inline Matrix operator^(const MatView& a, const MatView& b) { return impl::multNT(a,b); }

inline ZMatrix operator*(const ZMatrix& a,  const ZMatrix& b)  { return impl::multNN(a,b); }
inline ZMatrix operator*(const ZMatrix& a,  const ZMatView& b) { return impl::multNN(a,b); }
inline ZMatrix operator*(const ZMatView& a, const ZMatrix& b)  { return impl::multNN(a,b); }
inline ZMatrix operator*(const ZMatView& a, const ZMatView& b) { return impl::multNN(a,b); }
inline ZMatrix operator%(const ZMatrix& a,  const ZMatrix& b)  { return impl::multTN(a,b); }
inline ZMatrix operator%(const ZMatrix& a,  const ZMatView& b) { return impl::multTN(a,b); }
inline ZMatrix operator%(const ZMatView& a, const ZMatrix& b)  { return impl::multTN(a,b); }
inline ZMatrix operator%(const ZMatView& a, const ZMatView& b) { return impl::multTN(a,b); }
inline ZMatrix operator^(const ZMatrix& a,  const ZMatrix& b)  { return impl::multNT(a,b); }
inline ZMatrix operator^(const ZMatrix& a,  const ZMatView& b) { return impl::multNT(a,b); }
inline ZMatrix operator^(const ZMatView& a, const ZMatrix& b)  { return impl::multNT(a,b); }
inline ZMatrix operator^(const ZMatView& a, const ZMatView& b) { return impl::multNT(a,b); }

// operator*=
inline Matrix& operator*=(Matrix& a,  const Matrix& b)  { a = a*b; return a; }
inline Matrix& operator*=(Matrix& a,  const MatView& b) { a = a*b; return a; }

inline ZMatrix& operator*=(ZMatrix& a,  const ZMatrix& b)  { a = a*b; return a; }
inline ZMatrix& operator*=(ZMatrix& a,  const ZMatView& b) { a = a*b; return a; }

// operator % between vectors (which will return dot products).
inline double operator%(const VectorB& a, const VectorB& b) { return btas::dotc(a, b); }
inline double operator%(const VectorB& a, const VecView& b) { return btas::dotc(a, b); }
inline double operator%(const VecView& a, const VectorB& b) { return btas::dotc(a, b); }
inline double operator%(const VecView& a, const VecView& b) { return btas::dotc(a, b); }
inline std::complex<double> operator%(const ZVectorB& a, const ZVectorB& b) { return btas::dotc(a, b); }
inline std::complex<double> operator%(const ZVectorB& a, const ZVecView& b) { return btas::dotc(a, b); }
inline std::complex<double> operator%(const ZVecView& a, const ZVectorB& b) { return btas::dotc(a, b); }
inline std::complex<double> operator%(const ZVecView& a, const ZVecView& b) { return btas::dotc(a, b); }

// operator * and % between Matrix and VectorB
inline VectorB operator*(const Matrix& a, const VectorB& b)  { VectorB out(a.extent(0)); btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0}); return out; }
inline VectorB operator*(const Matrix& a, const VecView& b)  { VectorB out(a.extent(0)); btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0}); return out; }
inline VectorB operator*(const MatView& a, const VectorB& b) { VectorB out(a.extent(0)); btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0}); return out; }
inline VectorB operator*(const MatView& a, const VecView& b) { VectorB out(a.extent(0)); btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0}); return out; }
inline VectorB operator%(const Matrix& a, const VectorB& b)  { VectorB out(a.extent(1)); btas::contract(1.0, a, {1,0}, b, {1}, 0.0, out, {0}); return out; }
inline VectorB operator%(const Matrix& a, const VecView& b)  { VectorB out(a.extent(1)); btas::contract(1.0, a, {1,0}, b, {1}, 0.0, out, {0}); return out; }
inline VectorB operator%(const MatView& a, const VectorB& b) { VectorB out(a.extent(1)); btas::contract(1.0, a, {1,0}, b, {1}, 0.0, out, {0}); return out; }
inline VectorB operator%(const MatView& a, const VecView& b) { VectorB out(a.extent(1)); btas::contract(1.0, a, {1,0}, b, {1}, 0.0, out, {0}); return out; }
inline VectorB operator%(const VectorB& a, const Matrix& b)  { VectorB out(a.extent(1)); btas::contract(1.0, b, {1,0}, a, {1}, 0.0, out, {0}); return out; }
inline VectorB operator%(const VecView& a, const Matrix& b)  { VectorB out(a.extent(1)); btas::contract(1.0, b, {1,0}, a, {1}, 0.0, out, {0}); return out; }
inline VectorB operator%(const VectorB& a, const MatView& b) { VectorB out(a.extent(1)); btas::contract(1.0, b, {1,0}, a, {1}, 0.0, out, {0}); return out; }
inline VectorB operator%(const VecView& a, const MatView& b) { VectorB out(a.extent(1)); btas::contract(1.0, b, {1,0}, a, {1}, 0.0, out, {0}); return out; }
inline Matrix  operator^(const VectorB& a, const VectorB& b) { Matrix out(a.extent(0), b.extent(0)); /*btas::contract(1.0, a, {0}, b, {1}, 0.0, out, {0,1});*/ dger_(a.extent(0), b.extent(0), 1.0, a.data(), 1, b.data(), 1, out.data(), a.extent(0)); return out; }
inline Matrix  operator^(const VectorB& a, const VecView& b) { Matrix out(a.extent(0), b.extent(0)); /*btas::contract(1.0, a, {0}, b, {1}, 0.0, out, {0,1});*/ dger_(a.extent(0), b.extent(0), 1.0, a.data(), 1, b.data(), 1, out.data(), a.extent(0)); return out; }
inline Matrix  operator^(const VecView& a, const VectorB& b) { Matrix out(a.extent(0), b.extent(0)); /*btas::contract(1.0, a, {0}, b, {1}, 0.0, out, {0,1});*/ dger_(a.extent(0), b.extent(0), 1.0, a.data(), 1, b.data(), 1, out.data(), a.extent(0)); return out; }
inline Matrix  operator^(const VecView& a, const VecView& b) { Matrix out(a.extent(0), b.extent(0)); /*btas::contract(1.0, a, {0}, b, {1}, 0.0, out, {0,1});*/ dger_(a.extent(0), b.extent(0), 1.0, a.data(), 1, b.data(), 1, out.data(), a.extent(0)); return out; }

inline ZVectorB operator*(const ZMatrix& a, const ZVectorB& b)  { ZVectorB out(a.extent(0)); btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0}); return out; }
inline ZVectorB operator*(const ZMatrix& a, const ZVecView& b)  { ZVectorB out(a.extent(0)); btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0}); return out; }
inline ZVectorB operator*(const ZMatView& a, const ZVectorB& b) { ZVectorB out(a.extent(0)); btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0}); return out; }
inline ZVectorB operator*(const ZMatView& a, const ZVecView& b) { ZVectorB out(a.extent(0)); btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0}); return out; }
// TODO % and ^ operators require specification of complex conjugate

// operator* with scalar 
inline Matrix&  operator*=(Matrix& a, const double b)  { blas::scale_n(b, a.data(), a.size()); return a; }
inline MatView& operator*=(MatView& a, const double b) { blas::scale_n(b, a.data(), a.size()); return a; }
inline Matrix operator*(const Matrix& a, const double b)  { Matrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline Matrix operator*(const MatView& a, const double b) { Matrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline Matrix operator*(const double b, const Matrix& a)  { Matrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline Matrix operator*(const double b, const MatView& a) { Matrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline Matrix&  operator/=(Matrix& a, const double b)  { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline MatView& operator/=(MatView& a, const double b) { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline Matrix operator/(const Matrix& a, const double b)  { Matrix c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }
inline Matrix operator/(const MatView& a, const double b) { Matrix c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }

inline ZMatrix&  operator*=(ZMatrix& a, const double b)  { blas::scale_n(b, a.data(), a.size()); return a; }
inline ZMatView& operator*=(ZMatView& a, const double b) { blas::scale_n(b, a.data(), a.size()); return a; }
inline ZMatrix operator*(const ZMatrix& a, const double b)  { ZMatrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZMatrix operator*(const ZMatView& a, const double b) { ZMatrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZMatrix operator*(const double b, const ZMatrix& a)  { ZMatrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZMatrix operator*(const double b, const ZMatView& a) { ZMatrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZMatrix&  operator/=(ZMatrix& a, const double b)  { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline ZMatView& operator/=(ZMatView& a, const double b) { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline ZMatrix operator/(const ZMatrix& a, const double b)  { ZMatrix c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }
inline ZMatrix operator/(const ZMatView& a, const double b) { ZMatrix c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }
inline ZMatrix&  operator*=(ZMatrix& a, const std::complex<double> b)  { blas::scale_n(b, a.data(), a.size()); return a; }
inline ZMatView& operator*=(ZMatView& a, const std::complex<double> b) { blas::scale_n(b, a.data(), a.size()); return a; }
inline ZMatrix operator*(const ZMatrix& a, const std::complex<double> b)  { ZMatrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZMatrix operator*(const ZMatView& a, const std::complex<double> b) { ZMatrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZMatrix operator*(const std::complex<double> b, const ZMatrix& a)  { ZMatrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZMatrix operator*(const std::complex<double> b, const ZMatView& a) { ZMatrix c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZMatrix&  operator/=(ZMatrix& a, const std::complex<double> b)  { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline ZMatView& operator/=(ZMatView& a, const std::complex<double> b) { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline ZMatrix operator/(const ZMatrix& a, const std::complex<double> b)  { ZMatrix c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }
inline ZMatrix operator/(const ZMatView& a, const std::complex<double> b) { ZMatrix c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }

inline VectorB&  operator*=(VectorB& a, const double b)  { blas::scale_n(b, a.data(), a.size()); return a; }
inline VecView& operator*=(VecView& a, const double b) { blas::scale_n(b, a.data(), a.size()); return a; }
inline VectorB operator*(const VectorB& a, const double b)  { VectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline VectorB operator*(const VecView& a, const double b) { VectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline VectorB operator*(const double b, const VectorB& a)  { VectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline VectorB operator*(const double b, const VecView& a) { VectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline VectorB&  operator/=(VectorB& a, const double b)  { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline VecView& operator/=(VecView& a, const double b) { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline VectorB operator/(const VectorB& a, const double b)  { VectorB c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }
inline VectorB operator/(const VecView& a, const double b) { VectorB c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }

inline ZVectorB&  operator*=(ZVectorB& a, const double b)  { blas::scale_n(b, a.data(), a.size()); return a; }
inline ZVecView& operator*=(ZVecView& a, const double b) { blas::scale_n(b, a.data(), a.size()); return a; }
inline ZVectorB operator*(const ZVectorB& a, const double b)  { ZVectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZVectorB operator*(const ZVecView& a, const double b) { ZVectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZVectorB operator*(const double b, const ZVectorB& a)  { ZVectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZVectorB operator*(const double b, const ZVecView& a) { ZVectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZVectorB&  operator/=(ZVectorB& a, const double b)  { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline ZVecView& operator/=(ZVecView& a, const double b) { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline ZVectorB operator/(const ZVectorB& a, const double b)  { ZVectorB c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }
inline ZVectorB operator/(const ZVecView& a, const double b) { ZVectorB c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }
inline ZVectorB&  operator*=(ZVectorB& a, const std::complex<double> b)  { blas::scale_n(b, a.data(), a.size()); return a; }
inline ZVecView& operator*=(ZVecView& a, const std::complex<double> b) { blas::scale_n(b, a.data(), a.size()); return a; }
inline ZVectorB operator*(const ZVectorB& a, const std::complex<double> b)  { ZVectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZVectorB operator*(const ZVecView& a, const std::complex<double> b) { ZVectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZVectorB operator*(const std::complex<double> b, const ZVectorB& a)  { ZVectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZVectorB operator*(const std::complex<double> b, const ZVecView& a) { ZVectorB c(a); blas::scale_n(b, c.data(), c.size()); return c; }
inline ZVectorB&  operator/=(ZVectorB& a, const std::complex<double> b)  { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline ZVecView& operator/=(ZVecView& a, const std::complex<double> b) { blas::scale_n(1.0/b, a.data(), a.size()); return a; }
inline ZVectorB operator/(const ZVectorB& a, const std::complex<double> b)  { ZVectorB c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }
inline ZVectorB operator/(const ZVecView& a, const std::complex<double> b) { ZVectorB c(a); blas::scale_n(1.0/b, c.data(), c.size()); return c; }

}

#endif

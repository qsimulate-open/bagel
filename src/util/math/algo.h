//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: algo.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_MATH_ALGO_H
#define __SRC_MATH_ALGO_H

#include <array>
#include <complex>
#include <type_traits>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <src/util/f77.h>
#include <src/util/math/zquatev/zquatev.h>

namespace bagel {

extern void dcsrmm_(const char *transa, const int m, const int n, const int k, const double alpha, const double* adata,
                    const int* acols, const int* arind, const double* b, const int ldb, const double beta,
                    double* c, const int ldc);
extern void vdmul_(const int n, const double*, const double*, double*);

template <typename... Args>
auto zquatev(Args&&... args) -> decltype(ts::zquatev(std::forward<Args>(args)...)) {
  return ts::zquatev(std::forward<Args>(args)...);
}

namespace detail {
namespace {
  // conj function
  template<typename T> T conj(const T& a) { throw std::logic_error("detail::conj"); }
  template<> double conj(const double& a) { return a; }
  template<> std::complex<double> conj(const std::complex<double>& a) { return std::conj(a); }

  // real function
  template<typename T> double real(const T& a) { throw std::logic_error("detail::real"); }
  template<> double real(const double& a) { return a; }
  template<> double real(const std::complex<double>& a) { return a.real(); }

  // real function
  template<bool S, typename T, typename U>
  struct make_complex_impl { };
  template<typename T, typename U>
  struct make_complex_impl<true,T,U> { U operator()(const T& a, const T& b, U& out) { return out = U{a, b}; } };
  template<typename T, typename U>
  struct make_complex_impl<false,T,U> { U operator()(const T& a, const T& b, U& out) { throw std::logic_error("make_complex_impl"); return U(); } };
  template<typename T, typename U>
  static void make_complex(const T& a, const T& b, U& out) {
    make_complex_impl<std::is_same<U, std::complex<T>>::value,T,U>()(a, b, out);
  }

  // taylor expansion
  template<int M, int N>
  struct taylor {
    double operator()(const double& x, const std::array<double,N>& a) { return a[M-1] + x * taylor<M-1, N>()(x, a); }
  };
  template<int N>
  struct taylor<1, N> {
    double operator()(const double& x, std::array<double,N> a) { return a[0]; }
  };
  template<int N>
  double taylor_expansion(const double& x, const std::array<double,N>& a) {
    static_assert(N > 0, "illegal call of taylor_expansion");
    return taylor<N, N>()(x, a);
  }

} }

namespace blas {

// Transpose
template<typename T, typename U = T,
         class = typename std::enable_if< (std::is_same<double,T>::value || std::is_same<std::complex<double>,T>::value)
                                          and  std::is_convertible<U, T>::value >::type
        >
void transpose(const T* a, const int b, const int c, T* d, const U fac = static_cast<U>(1.0));
template<>
void transpose(const double* a, const int b, const int c, double* d, const double fac);
template<>
void transpose(const std::complex<double>* a, const int b, const int c, std::complex<double>* d, const std::complex<double> fac);
template<>
void transpose(const std::complex<double>* a, const int b, const int c, std::complex<double>* d, const double fac);

template<typename T, typename U = T,
         class = typename std::enable_if< (std::is_same<double,T>::value || std::is_same<std::complex<double>,T>::value)
                                          and  std::is_convertible<U, T>::value >::type
        >
void transpose_add(const T* a, const int b, const int c, T* d, const U fac = static_cast<U>(1.0));
template<>
void transpose_add(const double* a, const int b, const int c, double* d, const double fac);
template<>
void transpose_add(const std::complex<double>* a, const int b, const int c, std::complex<double>* d, const std::complex<double> fac);
template<>
void transpose_add(const std::complex<double>* a, const int b, const int c, std::complex<double>* d, const double fac);

template<typename T, typename U = T,
         class = typename std::enable_if< std::is_same<std::complex<double>,T>::value and std::is_convertible<U, T>::value >::type
        >
void transpose_conjg(const T* a, const int b, const int c, T* d, const U fac = static_cast<U>(1.0));
template<>
void transpose_conjg(const std::complex<double>* a, const int b, const int c, std::complex<double>* d, const std::complex<double> fac);
template<>
void transpose_conjg(const std::complex<double>* a, const int b, const int c, std::complex<double>* d, const double fac);


namespace {
  // AXPY
  template<class T, class U, typename Type,
           // T and U have to be either raw pointers or random access iterators
           class = typename std::enable_if< (std::is_pointer<T>::value) &&
                                            (std::is_pointer<U>::value) >::type >
  void ax_plus_y_n(const Type& a, const T p, const size_t n, U q) {
    std::transform(p, p+n, q, q, [&a](decltype(*p) i, decltype(*q) j) { return j+a*i; });
  }
  template<>
  void ax_plus_y_n(const double& a, const double* p, const size_t n, double* q) {
    daxpy_(n, a, p, 1, q, 1);
  }
  template<>
  void ax_plus_y_n(const std::complex<double>& a, const std::complex<double>* p, const size_t n, std::complex<double>* q) {
    zaxpy_(n, a, p, 1, q, 1);
  }
  template<>
  void ax_plus_y_n(const double& a, const std::complex<double>* p, const size_t n, std::complex<double>* q) {
    zaxpy_(n, static_cast<std::complex<double>>(a), p, 1, q, 1);
  }

  // DOT
  template<class T, class U,
           // T and U have to be either raw pointers or random access iterators
           class = typename std::enable_if< (std::is_pointer<T>::value) &&
                                            (std::is_pointer<U>::value) >::type >
  auto dot_product(const T p, const size_t n, const U q) -> decltype(*p * *q) {
    using ResultType = decltype(*p * *q);
    return std::inner_product(p, p+n, q, static_cast<ResultType>(0.0), std::plus<ResultType>(), [](decltype(*p) i, decltype(*q) j) { return detail::conj(i)*j; });
  }
  template<>
  double dot_product(const double* p, const size_t n, const double* q) {
    return ddot_(n, p, 1, q, 1);
  }
  template<>
  std::complex<double> dot_product(const std::complex<double>* p, const size_t n, const std::complex<double>* q) {
    return zdotc_(n, p, 1, q, 1);
  }

  template<class T, class U,
           // T and U have to be either raw pointers or random access iterators
           class = typename std::enable_if< (std::is_pointer<T>::value) &&
                                            (std::is_pointer<U>::value) >::type >
  auto dot_product_noconj(const T p, const size_t n, const U q) -> decltype(*p * *q) {
    using ResultType = decltype(*p * *q);
    return std::inner_product(p, p+n, q, static_cast<ResultType>(0.0), std::plus<ResultType>(), [](decltype(*p) i, decltype(*q) j) { return i*j; });
  }
  template<>
  double dot_product_noconj(const double* p, const size_t n, const double* q) {
    return ddot_(n, p, 1, q, 1);
  }
  template<>
  std::complex<double> dot_product_noconj(const std::complex<double>* p, const size_t n, const std::complex<double>* q) {
    return zdotu_(n, p, 1, q, 1);
  }


  // SCAL
  template<class T, class U,
           // U has to be either raw pointers or random access iterators
           class = typename std::enable_if<std::is_pointer<U>::value>::type >
  void scale_n(const T p, U q, const size_t n) {
    std::for_each(q, q+n, [&p](decltype(*q) j) { j *= p; });
  }
  template<>
  void scale_n(const double p, double* q, const size_t n) { dscal_(n, p, q, 1); }
  template<>
  void scale_n(const double p, std::complex<double>* q, const size_t n) { zscal_(n, p, q, 1); }
  template<>
  void scale_n(const std::complex<double> p, std::complex<double>* q, const size_t n) { zscal_(n, p, q, 1); }

  // average
  template<class T>
  auto average(const T& container) -> decltype(std::accumulate(container.begin(), container.end(), 0.0)/container.size()) {
    return std::accumulate(container.begin(), container.end(), 0.0) / container.size();
  }

  // conjugate an array
  template<class T,
           class = typename std::enable_if<std::is_pointer<T>::value>::type >
  void conj_n(T p, const size_t n) { throw std::logic_error("illegal call: blas::conj_n"); }
  template<>
  void conj_n(double* p, const size_t n) { /*do nothing*/ }
  template<>
  void conj_n(std::complex<double>* p, const size_t n) {
    double* dp = reinterpret_cast<double*>(p) + 1;
    for (double* i = dp; i <= dp + 2*n-2; i += 2) *i = -*i;
  }

}}

}

#endif

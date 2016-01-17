//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: matop.h
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

#ifndef __SRC_MATH_MATOP_H
#define __SRC_MATH_MATOP_H

#include <src/util/math/vectorb.h>
#include <src/util/math/matrix.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

namespace detail {
  template <typename T>
  struct is_mat {
    static const bool value = std::is_base_of<MatView, T>::value or std::is_base_of<Matrix, T>::value;
  };
  template <typename T>
  struct is_zmat {
    static const bool value = std::is_base_of<ZMatView, T>::value or std::is_base_of<ZMatrix, T>::value;
  };
  template <typename T>
  struct is_vec {
    static const bool value = std::is_base_of<VecView, T>::value or std::is_base_of<VectorB, T>::value;
  };
  template <typename T>
  struct is_zvec {
    static const bool value = std::is_base_of<ZVecView, T>::value or std::is_base_of<ZVectorB, T>::value;
  };
  template <typename T, typename U>
  struct is_matrix_pair {
    static const bool value = (is_mat<T>::value and is_mat<U>::value) or (is_zmat<T>::value and is_zmat<U>::value);
  };
  template <typename T, typename U>
  struct is_vector_pair {
    static const bool value = (is_vec<T>::value and is_vec<U>::value) or (is_zvec<T>::value and is_zvec<U>::value);
  };
  template <typename T, typename U>
  struct is_valid_pair {
    static const bool value = is_matrix_pair<T, U>::value or is_vector_pair<T, U>::value;
  };
  template <typename T>
  struct is_any_matrix {
    static const bool value = is_mat<T>::value or is_zmat<T>::value or is_vec<T>::value or is_zvec<T>::value;
  };
  template <typename T>
  struct returnable {
    using type = typename std::conditional<is_mat<T>::value, Matrix,
                   typename std::conditional<is_zmat<T>::value, ZMatrix,
                     typename std::conditional<is_vec<T>::value, VectorB,
                       typename std::conditional<is_zvec<T>::value, ZVectorB, void>::type>::type>::type>::type;
  };
}

// operator+= and -=
template <class T, class U,
          class = typename std::enable_if<detail::is_valid_pair<T,U>::value>::type
         >
T& operator+=(T& a,  const U& b) { assert(a.size() == b.size()); blas::ax_plus_y_n( 1.0, b.data(), a.size(), a.data()); return a; }

template <class T, class U,
          class = typename std::enable_if<detail::is_valid_pair<T,U>::value>::type
         >
T& operator-=(T& a,  const U& b) { assert(a.size() == b.size()); blas::ax_plus_y_n(-1.0, b.data(), a.size(), a.data()); return a; }

// operator+ and -
template <class T, class U,
          class = typename std::enable_if<detail::is_valid_pair<T,U>::value>::type
         >
typename detail::returnable<T>::type operator+(const T& a,  const U& b) { typename detail::returnable<T>::type out(a); out += b; return out; }

template <class T, class U,
          class = typename std::enable_if<detail::is_valid_pair<T,U>::value>::type
         >
typename detail::returnable<T>::type operator-(const T& a,  const U& b) { typename detail::returnable<T>::type out(a); out -= b; return out; }


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
    assert(a.mdim() == b.ndim());
    const int n = b.mdim();
    ZMatrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    const int m = a.mdim();
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      contract(1.0, a, {0,1}, b, {1,2}, 0.0, out, {0,2}, false, false);
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
    assert(a.ndim() == b.ndim());
    const int n = b.mdim();
    ZMatrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    const int m = a.ndim();
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      contract(1.0, a, {1,0}, b, {1,2}, 0.0, out, {0,2}, true, false);
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
    assert(a.mdim() == b.mdim());
    const int n = b.ndim();
    ZMatrix out(l, n, a.localized());

#ifdef HAVE_SCALAPACK
    const int m = a.mdim();
    assert(a.localized() == b.localized());
    if (a.localized() || std::min(std::min(l,m),n) < blocksize__) {
#endif
      contract(1.0, a, {0,1}, b, {2,1}, 0.0, out, {0,2}, false, true);
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
template <class T, class U,
          class = typename std::enable_if<detail::is_matrix_pair<T,U>::value>::type
         >
auto operator*(const T& a,  const U& b) ->decltype(impl::multNN(a,b)) { return impl::multNN(a,b); }

template <class T, class U,
          class = typename std::enable_if<detail::is_matrix_pair<T,U>::value>::type
         >
auto operator%(const T& a,  const U& b) ->decltype(impl::multTN(a,b)) { return impl::multTN(a,b); }

template <class T, class U,
          class = typename std::enable_if<detail::is_matrix_pair<T,U>::value>::type
         >
auto operator^(const T& a,  const U& b) ->decltype(impl::multNT(a,b)) { return impl::multNT(a,b); }

// operator*=
template <class T, class U,
          class = typename std::enable_if<detail::is_matrix_pair<T,U>::value and std::is_base_of<typename detail::returnable<T>::type, T>::value>::type
         >
T& operator*=(T& a,  const U& b) { a = a*b; return a; }

// operator % between vectors (which will return dot products).
template <typename T, typename U, class = typename std::enable_if<detail::is_vector_pair<T, U>::value>::type>
auto operator%(const T& a, const U& b) -> decltype(btas::dotc(a, b)) { return btas::dotc(a, b); }

// operator * and % between Matrix and VectorB
template <class T, class U,
          class = typename std::enable_if<(detail::is_mat<T>::value and detail::is_vec<U>::value) or (detail::is_zmat<T>::value and detail::is_zvec<U>::value)>::type
         >
typename detail::returnable<U>::type operator*(const T& a, const U& b)  {
  typename detail::returnable<U>::type out(a.extent(0));
  btas::contract(1.0, a, {0,1}, b, {1}, 0.0, out, {0});
  return out;
}

template <class T, class U,
          class = typename std::enable_if<(detail::is_mat<T>::value and detail::is_vec<U>::value) or (detail::is_zmat<T>::value and detail::is_zvec<U>::value)>::type
         >
typename detail::returnable<U>::type operator%(const T& a, const U& b)  {
  typename detail::returnable<U>::type out(a.extent(1));
  btas::contract(1.0, a, {1,0}, b, {1}, 0.0, out, {0}, true, false);
  return out;
}

template <class T, class U,
          class = typename std::enable_if<(detail::is_mat<T>::value and detail::is_vec<U>::value) or (detail::is_zmat<T>::value and detail::is_zvec<U>::value)>::type
         >
typename detail::returnable<U>::type operator%(const U& b, const T& a)  {
  auto out = a % b;
  conj_n(out.data(), out.size());
  return out;
}

// TODO to be cleaned up
inline Matrix  operator^(const VectorB& a, const VectorB& b) {
  Matrix out(a.extent(0), b.extent(0));
  /*btas::contract(1.0, a, {0}, b, {1}, 0.0, out, {0,1});*/
  dger_(a.extent(0), b.extent(0), 1.0, a.data(), 1, b.data(), 1, out.data(), a.extent(0));
  return out;
}
inline Matrix  operator^(const VectorB& a, const VecView& b) {
  Matrix out(a.extent(0), b.extent(0));
  /*btas::contract(1.0, a, {0}, b, {1}, 0.0, out, {0,1});*/
  dger_(a.extent(0), b.extent(0), 1.0, a.data(), 1, b.data(), 1, out.data(), a.extent(0));
  return out;
}
inline Matrix  operator^(const VecView& a, const VectorB& b) {
  Matrix out(a.extent(0), b.extent(0));
  /*btas::contract(1.0, a, {0}, b, {1}, 0.0, out, {0,1});*/
  dger_(a.extent(0), b.extent(0), 1.0, a.data(), 1, b.data(), 1, out.data(), a.extent(0));
  return out;
}
inline Matrix  operator^(const VecView& a, const VecView& b) {
  Matrix out(a.extent(0), b.extent(0));
  /*btas::contract(1.0, a, {0}, b, {1}, 0.0, out, {0,1});*/
  dger_(a.extent(0), b.extent(0), 1.0, a.data(), 1, b.data(), 1, out.data(), a.extent(0));
  return out;
}


// operator* with scalar
template <class T, typename U,
          class = typename std::enable_if<detail::is_any_matrix<T>::value and
                                          std::is_convertible<U, typename T::value_type>::value
                                         >::type
         >
T& operator*=(T& a, const U b) { blas::scale_n(b, a.data(), a.size()); return a; }

template <class T, typename U,
          class = typename std::enable_if<detail::is_any_matrix<T>::value and
                                          std::is_convertible<U, typename T::value_type>::value
                                         >::type
         >
T& operator/=(T& a, const U b) { a *= static_cast<typename T::value_type>(1.0)/b; return a; }

template <class T, typename U,
          class = typename std::enable_if<detail::is_any_matrix<T>::value and
                                          std::is_convertible<U, typename T::value_type>::value
                                         >::type
         >
typename detail::returnable<T>::type operator*(const T& a, const U b) { typename detail::returnable<T>::type c(a); c *= b; return c; }

template <class T, typename U,
          class = typename std::enable_if<detail::is_any_matrix<T>::value and
                                          std::is_convertible<U, typename T::value_type>::value
                                         >::type
         >
typename detail::returnable<T>::type operator*(const U b, const T& a) { return a*b; }

template <class T, typename U,
          class = typename std::enable_if<detail::is_any_matrix<T>::value and
                                          std::is_convertible<U, typename T::value_type>::value
                                         >::type
         >
typename detail::returnable<T>::type operator/(const T& a, const U b) { typename detail::returnable<T>::type c(a); c /= b; return c; }

}

#endif

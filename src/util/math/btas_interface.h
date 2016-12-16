//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: btas_interface.h
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

#ifndef __SRC_MATH_BTAS_INTERFACE_H
#define __SRC_MATH_BTAS_INTERFACE_H

#include <bagel_config.h>
#define _CONTRACT_OPT_BAGEL
#define _HAS_CBLAS
#ifdef HAVE_MKL_H
#define _HAS_INTEL_MKL
#endif

#include <iostream>
#include <iomanip>
#include <complex>
#include <btas/btas.h>
#include <btas/tensor.h>
#include <btas/tensor_func.h>
#include <src/util/math/btas_varray.h>
#include <src/util/math/preallocarray.h>

namespace btas {
  // int N is not nessesary, but leave it so that we can switch to fixed-rank tensors in the future
  template<int N>
  using CRange = RangeNd<CblasColMajor>;

  template<typename T>
  using Tensor1 = Tensor<T, CRange<1>, bagel::varray<T>>;
  template<typename T>
  using Tensor2 = Tensor<T, CRange<2>, bagel::varray<T>>;
  template<typename T>
  using Tensor3 = Tensor<T, CRange<3>, bagel::varray<T>>;
  template<typename T>
  using Tensor4 = Tensor<T, CRange<4>, bagel::varray<T>>;

  template<typename T>
  using TensorView1 = TensorView<T, CRange<1>, bagel::varray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
  template<typename T>
  using TensorView2 = TensorView<T, CRange<2>, bagel::varray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
  template<typename T>
  using TensorView3 = TensorView<T, CRange<3>, bagel::varray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
  template<typename T>
  using TensorView4 = TensorView<T, CRange<4>, bagel::varray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;

  template<typename T, int N>
  using TensorN = Tensor<T, CRange<N>, bagel::varray<T>>;
  template<typename T, int N>
  using TensorViewN = TensorView<T, CRange<N>, bagel::varray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;

  // Tensors on preallocated memory
  template<typename T>
  using PTensor1 = Tensor<T, CRange<1>, bagel::PreAllocArray<T>>;
  template<typename T>
  using PTensor2 = Tensor<T, CRange<2>, bagel::PreAllocArray<T>>;
  template<typename T>
  using PTensor3 = Tensor<T, CRange<3>, bagel::PreAllocArray<T>>;
  template<typename T>
  using PTensor4 = Tensor<T, CRange<4>, bagel::PreAllocArray<T>>;

  template<typename T>
  using PTensorView1 = TensorView<T, CRange<1>, bagel::PreAllocArray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
  template<typename T>
  using PTensorView2 = TensorView<T, CRange<2>, bagel::PreAllocArray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
  template<typename T>
  using PTensorView3 = TensorView<T, CRange<3>, bagel::PreAllocArray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
  template<typename T>
  using PTensorView4 = TensorView<T, CRange<4>, bagel::PreAllocArray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;

  template<typename T, int N>
  using PTensorN = Tensor<T, CRange<N>, bagel::PreAllocArray<T>>;
  template<typename T, int N>
  using PTensorViewN = TensorView<T, CRange<N>, bagel::PreAllocArray<T>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;

  // print functions for Tensor2
  template <typename T, class = typename std::enable_if<btas::is_boxtensor<T>::value>::type>
  void print(const T& t, std::string name = "", const typename T::size_type size = 10) {
    assert(t.rank() == 2 && t.range().ordinal().contiguous());
    if (!name.empty()) std::cout << "++++ " + name + " ++++" << std::endl;
    constexpr const int width = std::is_same<typename std::remove_cv<typename T::value_type>::type, double>::value ? 12 : 30;
    constexpr const int prec  = std::is_same<typename std::remove_cv<typename T::value_type>::type, double>::value ? 9 : 8;
    for (int i = 0; i != std::min(size, t.extent(0)); ++i) {
      for (int j = 0; j != std::min(size, t.extent(1)); ++j) {
        std::cout << std::fixed << std::setw(width) << std::setprecision(prec) << *(t.begin()+i+t.extent(0)*j)  << " ";
      }
      std::cout << std::endl;
    }
  }

}

extern template class btas::Tensor    <double,btas::RangeNd<CblasColMajor>,bagel::varray<double>>;
extern template class btas::TensorView<double,btas::RangeNd<CblasColMajor>,bagel::varray<double>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
extern template class btas::Tensor    <std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::varray<std::complex<double>>>;
extern template class btas::TensorView<std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::varray<std::complex<double>>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;

extern template class btas::Tensor    <double,btas::RangeNd<CblasColMajor>,bagel::PreAllocArray<double>>;
extern template class btas::TensorView<double,btas::RangeNd<CblasColMajor>,bagel::PreAllocArray<double>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
extern template class btas::Tensor    <std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::PreAllocArray<std::complex<double>>>;
extern template class btas::TensorView<std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::PreAllocArray<std::complex<double>>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(btas::Tensor1<double>)
BOOST_CLASS_EXPORT_KEY(btas::Tensor1<std::complex<double>>)

#endif

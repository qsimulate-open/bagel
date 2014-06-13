//
// BAGEL - Parallel electron correlation program.
// Filename: btas_interface.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_MATH_BTAS_INTERFACE_H
#define __SRC_MATH_BTAS_INTERFACE_H

#define _HAS_CBLAS
#define _HAS_INTEL_MKL

#include <btas/btas.h>
#include <btas/tensor.h>

namespace btas {
  template<typename T>
  using Tensor1 = Tensor<T, RangeNd<CblasColMajor, std::array<long,1>>, std::vector<T>>;
  template<typename T>
  using Tensor2 = Tensor<T, RangeNd<CblasColMajor, std::array<long,2>>, std::vector<T>>;
  template<typename T>
  using Tensor3 = Tensor<T, RangeNd<CblasColMajor, std::array<long,3>>, std::vector<T>>;
  template<typename T>
  using Tensor4 = Tensor<T, RangeNd<CblasColMajor, std::array<long,4>>, std::vector<T>>;

  template<typename T, int N>
  using TensorN = Tensor<T, RangeNd<CblasColMajor, std::array<long,N>>, std::vector<T>>;
}

#endif

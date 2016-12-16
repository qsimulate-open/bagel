//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: btas_interface.cc
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

#include <src/util/serialization.h>
#include <src/util/math/btas_interface.h>

template class btas::Tensor    <std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::varray<std::complex<double>>>;
template class btas::TensorView<std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::varray<std::complex<double>>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
template class btas::Tensor    <double,btas::RangeNd<CblasColMajor>,bagel::varray<double>>;
template class btas::TensorView<double,btas::RangeNd<CblasColMajor>,bagel::varray<double>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;

BOOST_CLASS_EXPORT_IMPLEMENT(btas::Tensor1<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(btas::Tensor1<std::complex<double>>)

template class btas::Tensor    <std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::PreAllocArray<std::complex<double>>>;
template class btas::TensorView<std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::PreAllocArray<std::complex<double>>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;
template class btas::Tensor    <double,btas::RangeNd<CblasColMajor>,bagel::PreAllocArray<double>>;
template class btas::TensorView<double,btas::RangeNd<CblasColMajor>,bagel::PreAllocArray<double>,btas::TensorViewPolicy<btas::TensorViewPolicy_RuntimeConst>>;

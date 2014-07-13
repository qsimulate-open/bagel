//
// BAGEL - Parallel electron correlation program.
// Filename: btas_interface.cc
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

#include <src/math/btas_interface.h>

template class btas::Tensor    <std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::varray<std::complex<double>>>;
template class btas::TensorView<std::complex<double>,btas::RangeNd<CblasColMajor>,bagel::varray<std::complex<double>>>;
template class btas::Tensor    <double,btas::RangeNd<CblasColMajor>,bagel::varray<double>>;
template class btas::TensorView<double,btas::RangeNd<CblasColMajor>,bagel::varray<double>>;

BOOST_CLASS_EXPORT_IMPLEMENT(btas::Tensor1<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(btas::Tensor1<std::complex<double>>)

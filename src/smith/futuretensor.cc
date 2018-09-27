//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: futuretensor.cc
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/futuretensor.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

template<typename DataType>
FutureTensor_<DataType>::FutureTensor_(const Tensor_<DataType>& i, shared_ptr<Task> j) : Tensor_<DataType>(i), init_(j) {
}

template<typename DataType>
void FutureTensor_<DataType>::init() const {
  init_->compute();
  initialized_ = true;
  mpi__->barrier();
}

template class bagel::SMITH::FutureTensor_<double>;
template class bagel::SMITH::FutureTensor_<complex<double>>;

#endif

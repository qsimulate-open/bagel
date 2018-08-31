//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: orthogonal_shift.cc
// Copyright (C) 2018 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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

#include <numeric>
#include <src/smith/moint.h>
#include <src/smith/orthogonal_basis.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


template<typename DataType>
void Orthogonal_Basis<DataType>::add_shift(shared_ptr<const Orthogonal_Basis<DataType>> t, const double shift, const bool imag) {
}


template<typename DataType>
tuple<shared_ptr<Matrix>,shared_ptr<Vec<double>>,shared_ptr<VecRDM<1>>,shared_ptr<VecRDM<2>>,shared_ptr<VecRDM<3>>,shared_ptr<VecRDM<3>>,vector<double>>
Orthogonal_Basis<DataType>::make_d2_imag(shared_ptr<const Orthogonal_Basis<DataType>> lambda, const double shift, const bool imag) const {
  auto dshift = make_shared<Matrix>(norb_, norb_);
  auto e0 = make_shared<Vec<double>>();
  auto e1 = make_shared<VecRDM<1>>();
  auto e2 = make_shared<VecRDM<2>>();
  auto e3 = make_shared<VecRDM<3>>();
  auto e4 = make_shared<VecRDM<3>>();
  vector<double> nimag;

  return tie(dshift, e0, e1, e2, e3, e4, nimag);
}

#endif

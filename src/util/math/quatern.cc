//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: quatern.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/util/math/quatern.h>

using namespace std;
using namespace bagel;

template<> void Quatern<double>::print() const {
  cout << setprecision(10) << setw(20) << data_[0];
  cout << setprecision(10) << setw(20) << data_[1];
  cout << setprecision(10) << setw(20) << data_[2];
  cout << setprecision(10) << setw(20) << data_[3];
  cout << endl;
};

template<> double Quatern<double>::norm() const { return sqrt(dot_product(*this)); };

template<> double Quatern<double>::dot_product(const Quatern<double>& o) const {
  return data_[0]*o.data_[0] + data_[1]*o.data_[1] + data_[2]*o.data_[2] + data_[3]*o.data_[3];
}

template<> void Quatern<double>::normalize() {
  const double n = norm();
  data_[0] /= n;
  data_[1] /= n;
  data_[2] /= n;
  data_[3] /= n;
}


//
// BAGEL - Parallel electron correlation program.
// Filename: quatern.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <src/math/quatern.h>

using namespace std;
using namespace bagel;

template<> void Quatern<double>::print() const {
  cout << setprecision(10) << setw(20) << data_[0];
  cout << setprecision(10) << setw(20) << data_[1];
  cout << setprecision(10) << setw(20) << data_[2];
  cout << setprecision(10) << setw(20) << data_[3];
  cout << endl;
};

template<> double Quatern<double>::norm() const { return ::sqrt(ddot(*this)); };

template<> double Quatern<double>::ddot(const Quatern<double>& o) const {
  return data_[0]*o.data_[0] + data_[1]*o.data_[1] + data_[2]*o.data_[2] + data_[3]*o.data_[3];
}

template<> void Quatern<double>::normalize() {
  const double n = norm();
  data_[0] /= n;
  data_[1] /= n;
  data_[2] /= n;
  data_[3] /= n;
};


//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: xyzfile.cc
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

#include <iomanip>
#include <iostream>
#include <src/util/math/xyzfile.h>

using namespace std;
using namespace bagel;


shared_ptr<XYZFile> XYZFile::clone() const {
  return make_shared<XYZFile>(mdim());
}


void XYZFile::print(const string in, const int dummy) const {
  cout << endl << "  * Nuclear energy gradient" << (in.empty() ? "" : (" " + in)) << endl << endl;
  for (int i = 0; i != mdim(); ++i) {
    cout << "    o Atom " << setw(3) << i << endl;
    cout << "        x  " << setprecision(10) << setw(20) << fixed << element(0,i) << endl;
    cout << "        y  " << setprecision(10) << setw(20) << fixed << element(1,i) << endl;
    cout << "        z  " << setprecision(10) << setw(20) << fixed << element(2,i) << endl;
  }
}

void XYZFile::print_export(const string in, const int dummy) const {
  auto matform = make_shared<Matrix>(mdim(), 3);
  for (int i = 0; i != mdim(); ++i) {
    matform->element(i, 0) = element(0, i);
    matform->element(i, 1) = element(1, i);
    matform->element(i, 2) = element(2, i);
  }
  matform->print();
}


shared_ptr<XYZFile> XYZFile::transform(const shared_ptr<const Matrix> a, const bool transpose) const {
  // a is ncart * ninternal quantity
  shared_ptr<XYZFile> out = clone();
  if (transpose) {
    dgemv_("T", a->ndim(), a->mdim(), 1.0, a->data(), a->ndim(), data(), 1, 0.0, out->data(), 1);
  } else {
    dgemv_("N", a->ndim(), a->mdim(), 1.0, a->data(), a->ndim(), data(), 1, 0.0, out->data(), 1);
  }
  return out;
}

double XYZFile::maximum(int atomno) const {
  double maximum = -1000.0;
  for (int i = 0; i != atomno; ++i) {
    if (fabs(element(0,i)) > maximum) maximum = fabs(element(0,i));
    if (fabs(element(1,i)) > maximum) maximum = fabs(element(1,i));
    if (fabs(element(2,i)) > maximum) maximum = fabs(element(2,i));
  }

  return maximum;
}

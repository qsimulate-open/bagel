//
// BAGEL - Parallel electron correlation program.
// Filename: xyzfile.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <iomanip>
#include <iostream>
#include <src/math/xyzfile.h>

using namespace std;
using namespace bagel;


shared_ptr<XYZFile> XYZFile::clone() const {
  return make_shared<XYZFile>(mdim_);
}


void XYZFile::print(const string in, const size_t dummy) const {
  cout << endl << "  * Nuclear energy gradient" << (in.empty() ? "" : (" " + in)) << endl << endl;
  for (int i = 0; i != mdim_; ++i) {
    cout << "    o Atom " << setw(3) << i << endl;
    cout << "        x  " << setprecision(10) << setw(20) << fixed << element(0,i) << endl;
    cout << "        y  " << setprecision(10) << setw(20) << fixed << element(1,i) << endl;
    cout << "        z  " << setprecision(10) << setw(20) << fixed << element(2,i) << endl;
  }
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

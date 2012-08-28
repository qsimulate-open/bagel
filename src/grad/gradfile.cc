//
// BAGEL - Parallel electron correlation program.
// Filename: gradfile.cc
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


#include <src/grad/gradfile.h>

using namespace std;

GradFile GradFile::operator+(const GradFile& o) const {
  GradFile out(*this);
  out.daxpy(1.0, o);
  return out;
}

GradFile GradFile::operator-(const GradFile& o) const {
  GradFile out(*this);
  out.daxpy(-1.0, o);
  return out;
}


shared_ptr<GradFile> GradFile::clone() const {
  const size_t natom = data_.size()/3;
  assert(data_.size()%3 == 0);
  shared_ptr<GradFile> out(new GradFile(natom));
  return out;
}

void GradFile::print() const {
  cout << endl << "  * Nuclear energy gradient" << endl << endl;
  assert(data_.size()%3 == 0);
  for (int i = 0; i != data_.size()/3; ++i) {
    cout << "    o Atom " << setw(3) << i << endl;
    cout << "        x  " << setprecision(10) << setw(20) << fixed << data(i,0) << endl;
    cout << "        y  " << setprecision(10) << setw(20) << fixed << data(i,1) << endl;
    cout << "        z  " << setprecision(10) << setw(20) << fixed << data(i,2) << endl;
  }
}

shared_ptr<GradFile> GradFile::transform(const unique_ptr<double[]>& a, const bool transpose) const {
  shared_ptr<GradFile> out = clone();
  const int size = data_.size();
  if (transpose) {
    dgemv_("T", size, size, 1.0, a.get(), size, data(), 1, 0.0, out->data(), 1);
  } else {
    dgemv_("N", size, size, 1.0, a.get(), size, data(), 1, 0.0, out->data(), 1);
  }
  return out;
}

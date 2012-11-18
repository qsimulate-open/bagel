//
// BAGEL - Parallel electron correlation program.
// Filename: rotfile.cc
// Copyright (C) 2011 Toru Shiozaki
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


#include <src/casscf/rotfile.h>

using namespace std;
using namespace bagel;

shared_ptr<RotFile> RotFile::clone() const {
  shared_ptr<RotFile> out(new RotFile(nclosed_, nact_, nvirt_, superci_));
  out->zero();
  return out;
}


shared_ptr<RotFile> RotFile::copy() const {
  shared_ptr<RotFile> out(new RotFile(*this));
  return out;
}


RotFile RotFile::operator+(const RotFile& o) const {
  RotFile out(*this);
  out.daxpy(1.0, o);
  return out;
}

RotFile RotFile::operator-(const RotFile& o) const {
  RotFile out(*this);
  out.daxpy(-1.0, o);
  return out;
}

RotFile& RotFile::operator+=(const RotFile& o) { daxpy(1.0, o); return *this; }
RotFile& RotFile::operator-=(const RotFile& o) { daxpy(-1.0, o); return *this; }

shared_ptr<Matrix> RotFile::unpack(shared_ptr<const Geometry> geom, const double a) const {

  const int nocc_ = nclosed_ + nact_;
  const int nbasis_ = nclosed_ + nact_ + nvirt_;
  shared_ptr<Matrix> out(new Matrix(nbasis_, nbasis_));
  fill(out->data(), out->data()+out->size(), a);

  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc_, i+nclosed_) = ele_va(j, i);
    }
    for (int j = 0; j != nclosed_; ++j) {
      out->element(i+nclosed_, j) = ele_ca(j, i);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc_, i) = ele_vc(j, i);
    }
  }
  for (int i = 0; i != nbasis_; ++i) {
    for (int j = 0; j <= i; ++j) {
      out->element(j, i) = -out->element(i, j);
    }
  }
  return out;
}


shared_ptr<Matrix> RotFile::unpack_sym(shared_ptr<const Geometry> geom, const double a) const {

  const int nocc_ = nclosed_ + nact_;
  const int nbasis_ = nclosed_ + nact_ + nvirt_;
  shared_ptr<Matrix> out(new Matrix(nbasis_, nbasis_));
  fill(out->data(), out->data()+out->size(), a);
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc_, i+nclosed_) = ele_va(j, i);
    }
    for (int j = 0; j != nclosed_; ++j) {
      out->element(i+nclosed_, j) = ele_ca(j, i);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc_, i) = ele_vc(j, i);
    }
  }
  for (int i = 0; i != nbasis_; ++i) {
    for (int j = 0; j <= i; ++j) {
      out->element(j, i) = out->element(i, j);
    }
  }
  return out;
}


double RotFile::orthog(list<shared_ptr<const RotFile> > c) {
  for (auto iter = c.begin(); iter != c.end(); ++iter)
    this->daxpy(- this->ddot(**iter), **iter);
  const double scal = 1.0/this->norm();
  dscal_(size_, scal, data(), 1);
  return 1.0/scal;
}

void RotFile::print() const {
  if (nact_ && nclosed_) {
    cout << " printing closed-active block" << endl;
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nclosed_; ++j) {
        cout << setw(10) << setprecision(6) << ele_ca(j,i);
      }
      cout << endl;
    }
  }
  if (nact_ && nvirt_) {
    cout << " printing virtual-active block" << endl;
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        cout << setw(10) << setprecision(6) << ele_va(j,i);
      }
      cout << endl;
    }
  }
  if (nclosed_ && nvirt_) {
    cout << " printing virtual-closed block" << endl;
    for (int i = 0; i != nclosed_; ++i) {
      for (int j = 0; j != nvirt_; ++j) {
        cout << setw(10) << setprecision(6) << ele_vc(j,i);
      }
      cout << endl;
    }
  }
  if (superci_) {
    cout << "reference weight " << ele_ref() << endl;
  }
}

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
// The BAGEL package is free software; you can redistribute it and/or modify
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


#include <iostream>
#include <iomanip>
#include <src/casscf/rotfile.h>

using namespace std;
using namespace bagel;


RotFile::RotFile(std::shared_ptr<const Matrix> o, const int iclos, const int iact, const int ivirt, const bool superci)
 : nclosed_(iclos), nact_(iact), nvirt_(ivirt), superci_(superci), size_(iclos*iact+iclos*ivirt+iact*ivirt+(superci ? 1 : 0)), data_(new double[size_]) {

  const int nocc = nclosed_ + nact_;
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      ele_va(j, i) = o->element(j+nocc, i+nclosed_);
    }
    for (int j = 0; j != nclosed_; ++j) {
      ele_ca(j, i) = o->element(i+nclosed_, j);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      ele_vc(j, i) = o->element(j+nocc, i);
    }
  }

}

shared_ptr<RotFile> RotFile::clone() const {
  auto out = make_shared<RotFile>(nclosed_, nact_, nvirt_, superci_);
  out->zero();
  return out;
}


shared_ptr<RotFile> RotFile::copy() const {
  return make_shared<RotFile>(*this);
}


shared_ptr<Matrix> RotFile::unpack(const double a) const {

  const int nocc = nclosed_ + nact_;
  const int nbasis = nclosed_ + nact_ + nvirt_;
  auto out = make_shared<Matrix>(nbasis, nbasis);
  fill_n(out->data(), out->size(), a);

  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc, i+nclosed_) = ele_va(j, i);
    }
    for (int j = 0; j != nclosed_; ++j) {
      out->element(i+nclosed_, j) = ele_ca(j, i);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc, i) = ele_vc(j, i);
    }
  }
  for (int i = 0; i != nbasis; ++i) {
    for (int j = 0; j <= i; ++j) {
      out->element(j, i) = -out->element(i, j);
    }
  }
  return out;
}


shared_ptr<Matrix> RotFile::unpack_sym(const double a) const {

  const int nocc = nclosed_ + nact_;
  const int nbasis = nclosed_ + nact_ + nvirt_;
  auto out = make_shared<Matrix>(nbasis, nbasis);
  fill_n(out->data(), out->size(), a);
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc, i+nclosed_) = ele_va(j, i);
    }
    for (int j = 0; j != nclosed_; ++j) {
      out->element(i+nclosed_, j) = ele_ca(j, i);
    }
  }
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_;   ++j) {
      out->element(j+nocc, i) = ele_vc(j, i);
    }
  }
  for (int i = 0; i != nbasis; ++i) {
    for (int j = 0; j <= i; ++j) {
      out->element(j, i) = out->element(i, j);
    }
  }
  return out;
}


double RotFile::orthog(list<shared_ptr<const RotFile>> c) {
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

//
// Newint - Parallel electron correlation program.
// Filename: rotfile.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/casscf/rotfile.h>

using namespace std;


shared_ptr<Matrix1e> RotFile::unpack(shared_ptr<Geometry> geom, const double a) const {

  const int nocc_ = nclosed_ + nact_;
  const int nbasis_ = nclosed_ + nact_ + nvirt_; 
  shared_ptr<Matrix1e> out(new Matrix1e(geom, nbasis_, nbasis_));
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

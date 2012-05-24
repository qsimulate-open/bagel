//
// Newint - Parallel electron correlation program.
// Filename: fock_base.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/fock.h>
#include <src/util/f77.h>
#include <src/scf/symmat.h>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <algorithm>

typedef std::shared_ptr<const Geometry> RefGeometry;
typedef std::shared_ptr<const Hcore> RefHcore;
typedef std::shared_ptr<Matrix1e> RefAODensity;
typedef std::shared_ptr<Shell> RefShell;
typedef std::shared_ptr<Fock_base> RefFock_base;

using namespace std;

Fock_base::Fock_base(const RefGeometry geom, const RefFock_base previous, const RefAODensity den, const vector<double>& schwarz)
 : Matrix1e(geom), previous_(previous), density_(den), schwarz_(schwarz) {

  schwarz_thresh_ = geom->schwarz_thresh();

  init(); // zero here
}


void Fock_base::fock_one_electron_part() {

  const int nirrep = geom_->nirrep();
  const int unit = 1;
  const double one = 1.0;
  const int size = nbasis_ * nbasis_;
  if (nirrep != 1) {
    Matrix1e intermediate(geom_);
    for (int i = 1; i != nirrep; ++i) {
      SymMat symm(geom_, i);
      Matrix1e tmp = symm % (*this) * symm;
      intermediate += tmp;
    }
    double* idata = intermediate.data();
    daxpy_(size, 1.0, idata, 1, data(), 1);
    dscal_(size, 1.0/nirrep, data(), 1);
  }
  const double* previous_data = previous_->data();
  daxpy_(size, 1.0, previous_data, 1, data(), 1); 

  fill_upper();
}


Fock_base::Fock_base(const RefGeometry geom, const RefHcore hcore)
 : Matrix1e(geom) {

  dcopy_(nbasis_*nbasis_, hcore->data(), 1, data(), 1); 

  fill_upper();
}


Fock_base::Fock_base(const RefGeometry geom)
 : Matrix1e(geom) {

  Hcore hcore(geom);
  dcopy_(nbasis_*nbasis_, hcore.data(), 1, data(), 1); 

  fill_upper();
}

Fock_base::~Fock_base() {

}


void Fock_base::computebatch(const vector<RefShell>& input, const int offsetb0, const int offsetb1, const int ndim) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb0 = input[1]->nbasis(); 
  const int dimb1 = input[0]->nbasis(); 

  for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j) {
      data_[i * ndim + j] = 0.0; 
    }
  }
}



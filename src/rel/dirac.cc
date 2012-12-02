//
// BAGEL - Parallel electron correlation program.
// Filename: dirac.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/util/constants.h>
#include <src/rel/dirac.h>
#include <src/util/zmatrix.h>
#include <src/util/matrix.h>

using namespace std;
using namespace bagel;

void Dirac::compute() {
  const int n = geom_->nbasis();
  const double c = c__;
  const double w = 0.25/(c*c);

  kinetic_->print("kinetic");
  nai_->print("nai");

  shared_ptr<Matrix> hcore(new Matrix(2*n, 2*n));

  // RKB hcore: T is off diagonal block matrices, V is first main diagonal, and 1/4m^2c^2W-T is second main diagonal
  hcore->copy_block(0, 0, n, n, nai_);
  hcore->copy_block(0, n, n, n, kinetic_);
  hcore->copy_block(n, 0, n, n, kinetic_);
  hcore->copy_block(n, n, n, n, shared_ptr<Matrix>(new Matrix(*nai_*w - *kinetic_)));
  
  unique_ptr<double[]> eig(new double[hcore->ndim()]);
  hcore->diagonalize(eig.get()); 
  
  hcore->print("hcore");

  for (int i = 0; i != 2*n; ++i) cout << setprecision(10) << setw(15) << eig[i] << endl;
}


shared_ptr<Reference> Dirac::conv_to_ref() const {
  assert(false);
  return shared_ptr<Reference>();
}

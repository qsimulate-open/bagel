//
// BAGEL - Parallel electron correlation program.
// Filename: fock_base.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/fock.h>
#include <src/scf/symmat.h>

using namespace std;
using namespace bagel;


Fock_base::Fock_base(const shared_ptr<const Geometry> geom, const shared_ptr<const Matrix> previous, const std::shared_ptr<const Matrix> den, const vector<double>& schwarz)
 : Matrix1e(geom), geom_(geom), previous_(previous), density_(den), schwarz_(schwarz) {

  schwarz_thresh_ = geom->schwarz_thresh();

  init(); // zero here
}


void Fock_base::fock_one_electron_part() {

  const int nbasis = ndim_;
  assert(ndim_ == mdim_);

  const int nirrep = geom_->nirrep();
  const int unit = 1;
  const double one = 1.0;
  const int size = nbasis * nbasis;
  if (nirrep != 1) {
    Matrix intermediate(nbasis, nbasis);
    for (int i = 1; i != nirrep; ++i) {
      SymMat symm(geom_, i);
      Matrix tmp = symm % (*this) * symm;
      intermediate += tmp;
    }
    *this += intermediate;
    *this *= 1.0/nirrep;
  }

  *this += *previous_;

  fill_upper();
}


void Fock_base::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb0 = input[1]->nbasis();
  const int dimb1 = input[0]->nbasis();

  for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j) {
      data_[i*ndim_+j] = 0.0;
    }
  }
}



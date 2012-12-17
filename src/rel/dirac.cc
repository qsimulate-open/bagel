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

  shared_ptr<ZMatrix> hcore = hcore_construct(n);
  shared_ptr<ZMatrix> s12 = s12_construct(n);

  ZMatrix interm = *s12 % *hcore * *s12;
  unique_ptr<double[]> eig(new double[hcore->ndim()]);
  interm.diagonalize(eig.get()); 

  print_eig(n, eig);
}

//Print non dirac sea eigenvalues
void Dirac::print_eig(const int n, const unique_ptr<double[]>& eig) {
  for (int i = 2*n; i != 4*n; ++i) cout << setprecision(10) << setw(15) << eig[i] << endl;
}

shared_ptr<Reference> Dirac::conv_to_ref() const {
  assert(false);
  return shared_ptr<Reference>();
}

shared_ptr<ZMatrix> Dirac::hcore_construct(const int n) {
  const double c = c__;
  const complex<double> w(0.25/(c*c), 0.0);

  const complex<double> coeff1 (1.0, 0.0);
  const complex<double> coeffi (0.0, 1.0);

  shared_ptr<ZMatrix> out(new ZMatrix(4*n, 4*n));
  shared_ptr<ZMatrix> znai(new ZMatrix(2*n, 2*n));
  shared_ptr<ZMatrix> zkinetic(new ZMatrix(2*n, 2*n));

  array<shared_ptr<ZMatrix>,4> zsmallnai;
  for (auto& i : zsmallnai)
    i = znai->clone(); 

  znai->copy_real_block(coeff1, 0, 0, n, n, nai_);
  znai->copy_real_block(coeff1, n, n, n, n, nai_);
  zkinetic->copy_real_block(coeff1, 0, 0, n, n, kinetic_);
  zkinetic->copy_real_block(coeff1, n, n, n, n, kinetic_);

  zsmallnai[0]->copy_real_block(coeff1, 0, 0, n, n, (*smallnai_)[0]);
  zsmallnai[0]->copy_real_block(coeff1, n, n, n, n, (*smallnai_)[0]);
  zsmallnai[1]->copy_real_block(coeffi, 0, 0, n, n, (*smallnai_)[1]);
  zsmallnai[1]->copy_real_block(-coeffi, n, n, n, n, (*smallnai_)[1]);
  zsmallnai[2]->copy_real_block(coeffi, 0, n, n, n, (*smallnai_)[2]);
  zsmallnai[2]->copy_real_block(coeffi, n, 0, n, n, (*smallnai_)[2]);
  zsmallnai[3]->copy_real_block(coeff1, 0, n, n, n, (*smallnai_)[3]);
  zsmallnai[3]->copy_real_block(-coeff1, n, 0, n, n, (*smallnai_)[3]);

  shared_ptr<ZMatrix> smallnai(new ZMatrix(*zsmallnai[0] + *zsmallnai[1] + *zsmallnai[2] + *zsmallnai[3]));

  // RKB hcore: T is off diagonal block matrices, V is first main diagonal, and 1/4m^2c^2W-T is second main diagonal
  out->zero();
  out->copy_block(0, 0, 2*n, 2*n, znai);
  out->copy_block(0, 2*n, 2*n, 2*n, zkinetic);
  out->copy_block(2*n, 0, 2*n, 2*n, zkinetic);
  out->copy_block(2*n, 2*n, 2*n, 2*n, shared_ptr<ZMatrix>(new ZMatrix(*smallnai * w - *zkinetic)));

  return out;
}

shared_ptr<ZMatrix> Dirac::s12_construct(const int n) {
  const double c = c__;
  const double t = 0.5/(c*c);
  const complex<double> coeff1 (1.0, 0.0);

  shared_ptr<ZMatrix> out(new ZMatrix(4*n, 4*n));
  shared_ptr<Overlap> ovl(new Overlap(*overlap_));

  shared_ptr<Matrix> k12(new Matrix(*kinetic_ * t));
  ovl->inverse_half();
  k12->inverse_half();

  out->zero();
  out->copy_real_block(coeff1, 0, 0, n, n, ovl);
  out->copy_real_block(coeff1, n, n, n, n, ovl);
  out->copy_real_block(coeff1, 2*n, 2*n, n, n, k12); 
  out->copy_real_block(coeff1, 3*n, 3*n, n, n, k12); 

  return out;
}

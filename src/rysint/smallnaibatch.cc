//
// BAGEL - Parallel electron correlation program.
// Filename: smallnaibatch.cc
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


#include <iostream>
#include <iomanip>
#include <src/rysint/naibatch.h>
#include <src/rysint/smallnaibatch.h>

using namespace std;
using namespace bagel;

SmallNAIBatch::SmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Geometry> geom)
  : geom_(geom), shells_(info), size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()) {

  for (int i = 0; i != 4; ++i)
     data_[i] = shared_ptr<Matrix>(new Matrix(shells_[0]->nbasis(), shells_[1]->nbasis(), true));
}


SmallNAIBatch::~SmallNAIBatch() {
}

void SmallNAIBatch::compute() {
  // first compute uncontracted NAI with auxiliary basis (cartesian)
  const shared_ptr<const Matrix> nai = nai_compute();

  std::array<shared_ptr<Matrix>,3> ints;
  for (int i = 0; i != 3; ++i)
    ints[i] = shared_ptr<Matrix>(new Matrix(*shells_[0]->small(i) % *nai));

  array<int,3> f = {{2,3,1}};
  array<int,3> b = {{3,1,2}};

  // 0) x^x + y^y + z^z
  // 1) x^y - y^x
  // 2) y^z - z^y
  // 3) z^x - x^z

  // -1 because <m|p|n>^dagger = -<n|p|m>  (can be proven by integration by part)
  for (int i = 0; i != 3; ++i) {
    *data_[0]    += *ints[i]      * *shells_[1]->small(i);
    *data_[b[i]] += *ints[b[i]-1] * *shells_[1]->small(i);
    *data_[i+1]  -= *ints[f[i]-1] * *shells_[1]->small(i);
  }

}


shared_ptr<Matrix> SmallNAIBatch::nai_compute() const {

  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int a0size_inc = shells_[0]->aux_inc()->nbasis();
  const int a1size_inc = shells_[1]->aux_inc()->nbasis();
  const int a0size_dec = shells_[0]->aux_dec() ? shells_[0]->aux_dec()->nbasis() : 0;
  const int a1size_dec = shells_[1]->aux_dec() ? shells_[1]->aux_dec()->nbasis() : 0;
  const int a0 = a0size_inc + a0size_dec;
  const int a1 = a1size_inc + a1size_dec;

  shared_ptr<Matrix> nai(new Matrix(a0, a1, true));
  {
    shared_ptr<NAIBatch> naic(new NAIBatch(array<shared_ptr<const Shell>,2>{{shells_[0]->aux_inc(), shells_[1]->aux_inc()}}, geom_));
    naic->compute();
    nai->copy_block(0, 0, a0size_inc, a1size_inc, naic->data());
  }
  if (shells_[0]->aux_dec() && shells_[1]->aux_dec()) {
    shared_ptr<NAIBatch> naic(new NAIBatch(array<shared_ptr<const Shell>,2>{{shells_[0]->aux_dec(), shells_[1]->aux_dec()}}, geom_));
    naic->compute();
    nai->copy_block(a0size_inc, a1size_inc, a0size_dec, a1size_dec, naic->data());
  }
  if (shells_[0]->aux_dec()) {
    shared_ptr<NAIBatch> naic(new NAIBatch(array<shared_ptr<const Shell>,2>{{shells_[0]->aux_dec(), shells_[1]->aux_inc()}}, geom_));
    naic->compute();
    nai->copy_block(a0size_inc, 0, a0size_dec, a1size_inc, naic->data());
  }
  if (shells_[1]->aux_dec()) {
    shared_ptr<NAIBatch> naic(new NAIBatch(array<shared_ptr<const Shell>,2>{{shells_[0]->aux_inc(), shells_[1]->aux_dec()}}, geom_));
    naic->compute();
    nai->copy_block(0, a1size_inc, a0size_inc, a1size_dec, naic->data());
  }
  return nai;
}

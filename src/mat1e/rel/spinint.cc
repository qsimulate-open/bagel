//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spinint.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <src/util/constants.h>
#include <src/mat1e/overlap.h>
#include <src/mat1e/rel/spinint.h>
#include <src/mat1e/rel/small1e_general.h>
#include <src/integral/os/overlapbatch.h>

using namespace std;
using namespace bagel;

void RelSpinInt::compute_() {
  const int n = geom_->nbasis();
  const complex<double> imag(0.0, 1.0);
  auto overlap = make_shared<Overlap>(geom_);

  for (int i = 0; i != 3; ++i) data_[i] = make_shared<ZMatrix>(4*n, 4*n);

  // Large component
  data_[0]->add_real_block( 0.5,   n,   0, n, n, *overlap);
  data_[0]->add_real_block( 0.5,   0,   n, n, n, *overlap);

  data_[1]->add_real_block( 0.5*imag,   n,   0, n, n, *overlap);
  data_[1]->add_real_block(-0.5*imag,   0,   n, n, n, *overlap);

  data_[2]->add_real_block( 0.5,   0,   0, n, n, *overlap);
  data_[2]->add_real_block(-0.5,   n,   n, n, n, *overlap);

  // Small component
  // TODO Simplify this code
  // Commented out lines cancel at zero-field; can be replaced with field * overlap for GIAO-RMB
  const double w = 1.0/(8.0*c__*c__);
  auto smallints = make_shared<Small1e_General<OverlapBatch>>(geom_);

  // x^x contributions
  data_[0]->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[0]);
  data_[0]->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[0]);
  data_[1]->add_real_block( imag*w, 2*n, 3*n, n, n, (*smallints)[0]);
  data_[1]->add_real_block(-imag*w, 3*n, 2*n, n, n, (*smallints)[0]);
  data_[2]->add_real_block(     -w, 2*n, 2*n, n, n, (*smallints)[0]);
  data_[2]->add_real_block(      w, 3*n, 3*n, n, n, (*smallints)[0]);

  // y^y contributions
  data_[0]->add_real_block(     -w, 2*n, 3*n, n, n, (*smallints)[1]);
  data_[0]->add_real_block(     -w, 3*n, 2*n, n, n, (*smallints)[1]);
  data_[1]->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[1]);
  data_[1]->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[1]);
  data_[2]->add_real_block(     -w, 2*n, 2*n, n, n, (*smallints)[1]);
  data_[2]->add_real_block(      w, 3*n, 3*n, n, n, (*smallints)[1]);

  // z^z contributions
  data_[0]->add_real_block(     -w, 2*n, 3*n, n, n, (*smallints)[2]);
  data_[0]->add_real_block(     -w, 3*n, 2*n, n, n, (*smallints)[2]);
  data_[1]->add_real_block( imag*w, 2*n, 3*n, n, n, (*smallints)[2]);
  data_[1]->add_real_block(-imag*w, 3*n, 2*n, n, n, (*smallints)[2]);
  data_[2]->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[2]);
  data_[2]->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[2]);

  // x^y contributions
  data_[0]->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[3]);
  data_[0]->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[3]);
  data_[1]->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[3]);
  data_[1]->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[3]);
  //data_[2]->add_real_block(-imag*w, 2*n, 2*n, n, n, (*smallints)[3]);
  //data_[2]->add_real_block(-imag*w, 3*n, 3*n, n, n, (*smallints)[3]);

  // y^x contributions
  data_[0]->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[6]);
  data_[0]->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[6]);
  data_[1]->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[6]);
  data_[1]->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[6]);
  //data_[2]->add_real_block( imag*w, 2*n, 2*n, n, n, (*smallints)[6]);
  //data_[2]->add_real_block( imag*w, 3*n, 3*n, n, n, (*smallints)[6]);

  // y^z contributions
  //data_[0]->add_real_block(-imag*w, 2*n, 2*n, n, n, (*smallints)[4]);
  //data_[0]->add_real_block(-imag*w, 3*n, 3*n, n, n, (*smallints)[4]);
  data_[1]->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[4]);
  data_[1]->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[4]);
  data_[2]->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[4]);
  data_[2]->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[4]);

  // z^y contributions
  //data_[0]->add_real_block( imag*w, 2*n, 2*n, n, n, (*smallints)[7]);
  //data_[0]->add_real_block( imag*w, 3*n, 3*n, n, n, (*smallints)[7]);
  data_[1]->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[7]);
  data_[1]->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[7]);
  data_[2]->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[7]);
  data_[2]->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[7]);

  // z^x contributions
  data_[0]->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[5]);
  data_[0]->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[5]);
  //data_[1]->add_real_block(-imag*w, 2*n, 2*n, n, n, (*smallints)[5]);
  //data_[1]->add_real_block(-imag*w, 3*n, 3*n, n, n, (*smallints)[5]);
  data_[2]->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[5]);
  data_[2]->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[5]);

  // x^z contributions
  data_[0]->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[8]);
  data_[0]->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[8]);
  //data_[1]->add_real_block( imag*w, 2*n, 2*n, n, n, (*smallints)[8]);
  //data_[1]->add_real_block( imag*w, 3*n, 3*n, n, n, (*smallints)[8]);
  data_[2]->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[8]);
  data_[2]->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[8]);
}


#ifndef NDEBUG
void RelTRevInt::compute_() {
  const int n = geom_->nbasis();
  copy_real_block(-1.0, 0, n, n, n, *overlap_);
  copy_real_block( 1.0, n, 0, n, n, *overlap_);
  copy_real_block(-0.5/(c__*c__), 2*n, 3*n, n, n, *kinetic_);
  copy_real_block( 0.5/(c__*c__), 3*n, 2*n, n, n, *kinetic_);
}
#endif



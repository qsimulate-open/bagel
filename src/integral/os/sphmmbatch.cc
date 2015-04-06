//
// BAGEL - Parallel electron correlation program.
// Filename: sphmmbatch.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#include <src/integral/carsphlist.h>
#include <src/integral/os/sphmmbatch.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;

SphMMBatch::SphMMBatch(const array<shared_ptr<const Shell>,2>& basis, shared_ptr<const Molecule> m, const int lmax)
 : basisinfo_(basis), mol_(m), lmax_(lmax)  {

  stack_ = resources__->get();

  size_block_ = basisinfo_[0]->nbasis() * basisinfo_[1]->nbasis();
  size_alloc_ = size_block_ * nblocks();
  stack_save_ = stack_->get(size_alloc_);
  data_ = stack_save_;

}


void SphMMBatch::compute() {

  const int dimb1 = basisinfo_[0]->nbasis();
  const int dimb0 = basisinfo_[1]->nbasis();

  const int ncblocks = lmax_*(lmax_*lmax_+6*lmax_+11)/6;
  vector<const double*> dat(ncblocks);


  MMBatch carmultipole(basisinfo_, mol_, lmax_);
  assert(ncblocks == carmultipole.num_blocks());
  carmultipole.compute();

  for (int i = 0; i != ncblocks; ++i)
    dat[i] = carmultipole.data() + carmultipole.size_block()*i;

  double car[100];
  double sph[100];
  assert(ncblocks<=100 && nblocks()<=100);
  for (int i = 0; i != dimb0; ++i) {
    for (int j = 0; j != dimb1; ++j) {
      double *tmp0, *tmp1;
      tmp0 = car;
      tmp1 = sph;

      for (int k = 0; k != ncblocks; ++k)
        car[k] = *dat[k]++;

      int offset0 = 0;
      int offset1 = 0;
      for (int l = 1; l <= lmax_; ++l) {
        const unsigned int carsph_index = l * ANG_HRR_END;
        const int nloops = 1;
        carsphlist.carsphfunc_call(carsph_index, nloops, tmp0, tmp1);
        offset0 += (l+1)*(l+2)/2;
        offset1 += 2*l+1;
        tmp0 = &car[offset0];
        tmp1 = &sph[offset1];
      }

      for (int k = 0; k != nblocks(); ++k) {
        data_[size_block_ * k + j + i * dimb1] = sph[k];
      }

    }
  }

}


SphMMBatch::~SphMMBatch() {
  stack_->release(size_alloc_, stack_save_);
  resources__->release(stack_);
}


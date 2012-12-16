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
#include <src/rel/smallnaibatch.h>
#include <src/osint/momentbatch.h>
#include <src/osint/overlapbatch.h>
#include <src/rysint/naibatch.h>
#include <src/rysint/carsphlist.h>
#include <src/util/matrix.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


SmallNAIBatch::SmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Geometry> geom)
  : geom_(geom), shells_(info), aux_inc_{{shells_[0]->kinetic_balance_uncont(1), shells_[1]->kinetic_balance_uncont(1)}},
    size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()), stack_(resources__->get()) {

  aux_dec_ = array<shared_ptr<const Shell>,2>{{shells_[0]->kinetic_balance_uncont(-1), shells_[1]->kinetic_balance_uncont(-1)}};

  for (int i = 0; i != 4; ++i)
     data_[i] = shared_ptr<Matrix>(new Matrix(shells_[0]->nbasis(), shells_[1]->nbasis()));
}


SmallNAIBatch::~SmallNAIBatch() {
  resources__->release(stack_);
}

void SmallNAIBatch::compute() {
  // then we need to have momentum integrals
  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int a0size_inc = aux_inc_[0]->nbasis();
  const int a1size_inc = aux_inc_[1]->nbasis();
  const int a0size_dec = aux_dec_[0] ? aux_dec_[0]->nbasis() : 0;
  const int a1size_dec = aux_dec_[1] ? aux_dec_[1]->nbasis() : 0;
  const int a0 = a0size_inc + a0size_dec;
  const int a1 = a1size_inc + a1size_dec;

  // first compute uncontracted NAI with auxiliary basis (cartesian)
#define LOCAL_DEBUG
#ifdef LOCAL_DEBUG
  const shared_ptr<const Matrix> nai = nai_compute();

#else
  assert(false);
#endif

#if 1
  std::array<shared_ptr<Matrix>,3> ints;

  shared_ptr<Matrix> ovl0 = ovl_compute(0);
  shared_ptr<Matrix> ovl1 = ovl_compute(1);

  {
    array<shared_ptr<Matrix>,3> tmp = moment_compute(0, ovl0);

    for (int i = 0; i != 3; ++i)
      ints[i] = shared_ptr<Matrix>(new Matrix(*tmp[i] % *nai));
  }

  array<int,3> f = {{2,3,1}};
  array<int,3> b = {{3,1,2}};

  {
    array<shared_ptr<Matrix>,3> tmp2 = moment_compute(1, ovl1);
    // 0) x^x + y^y + z^z
    // 1) x^y - y^x
    // 2) y^z - z^y
    // 3) z^x - x^z

    // -1 because <m|p|n>^dagger = -<n|p|m>  (can be proven by integration by part)
    for (int i = 0; i != 3; ++i) {
      *data_[0]    += *ints[i] * *tmp2[i];
      *data_[b[i]] += *ints[b[i]-1] * *tmp2[i];
      *data_[i+1]  -= *ints[f[i]-1] * *tmp2[i];
    }

  }

#endif
}

shared_ptr<Matrix> SmallNAIBatch::ovl_compute(const int num) const {

  const int asize_inc = aux_inc_[num]->nbasis();
  const int asize_dec = aux_dec_[num] ? aux_dec_[num]->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  std::shared_ptr<Matrix> overlap(new Matrix(a,a));

  {
    OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_inc_[num], aux_inc_[num]}}, stack_);
    ovl.compute();
    for (int i = 0; i != aux_inc_[num]->nbasis(); ++i)
      copy_n(ovl.data() + i*(aux_inc_[num]->nbasis()), aux_inc_[num]->nbasis(), overlap->data() + i*a);
  }
  if (aux_dec_[num]) {
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_dec_[num], aux_dec_[num]}}, stack_);
      ovl.compute();
      for (int i = aux_inc_[num]->nbasis(), j = 0; i != a; ++i, ++j) 
        copy_n(ovl.data() + j*(aux_dec_[num]->nbasis()), aux_dec_[num]->nbasis(), overlap->data() + i*a + aux_inc_[num]->nbasis());
    }
    {
      OverlapBatch ovl(array<shared_ptr<const Shell>,2>{{aux_inc_[num], aux_dec_[num]}}, stack_);
      ovl.compute();
      for (int i = aux_inc_[num]->nbasis(); i != a; ++i)
        for (int j = 0; j != aux_inc_[num]->nbasis(); ++j)
          overlap->element(j,i) = overlap->element(i,j) = *(ovl.data()+j+aux_inc_[num]->nbasis()*(i-aux_inc_[num]->nbasis()));
    }
  }
  overlap->inverse();
  return overlap;
}

shared_ptr<Matrix> SmallNAIBatch::nai_compute() const {

  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int a0size_inc = aux_inc_[0]->nbasis();
  const int a1size_inc = aux_inc_[1]->nbasis();
  const int a0size_dec = aux_dec_[0] ? aux_dec_[0]->nbasis() : 0;
  const int a1size_dec = aux_dec_[1] ? aux_dec_[1]->nbasis() : 0;
  const int a0 = a0size_inc + a0size_dec;
  const int a1 = a1size_inc + a1size_dec;

  shared_ptr<Matrix> nai(new Matrix(a0, a1));
  {
    shared_ptr<NAIBatch> naic(new NAIBatch(aux_inc_, geom_, stack_));
    naic->compute();
    nai->copy_block(0, 0, a0size_inc, a1size_inc, naic->data());
  }
  if (aux_dec_[0] && aux_dec_[1]) {
    shared_ptr<NAIBatch> naic(new NAIBatch(aux_dec_, geom_, stack_));
    naic->compute();
    nai->copy_block(a0size_inc, a1size_inc, a0size_dec, a1size_dec, naic->data());
  }
  if (aux_dec_[0]) {
    shared_ptr<NAIBatch> naic(new NAIBatch(array<shared_ptr<const Shell>,2>{{aux_dec_[0],aux_inc_[1]}}, geom_, stack_));
    naic->compute();
    nai->copy_block(a0size_inc, 0, a0size_dec, a1size_inc, naic->data());
  }
  if (aux_dec_[1]) {
    shared_ptr<NAIBatch> naic(new NAIBatch(array<shared_ptr<const Shell>,2>{{aux_inc_[0],aux_dec_[1]}}, geom_, stack_));
    naic->compute();
    nai->copy_block(0, a1size_inc, a0size_inc, a1size_dec, naic->data());
  }
  return nai;
}


array<shared_ptr<Matrix>,3>  SmallNAIBatch::moment_compute(const int num, const shared_ptr<const Matrix> overlap) const {
  const int ssize = shells_[num]->nbasis();
  const int asize_inc = aux_inc_[num]->nbasis();
  const int asize_dec = aux_dec_[num] ? aux_dec_[num]->nbasis() : 0;
  const int a = asize_inc + asize_dec;

  shared_ptr<MomentBatch> coeff0(new MomentBatch(array<shared_ptr<const Shell>,2>{{shells_[num]->cartesian_shell(), aux_inc_[num]}}, stack_));
  coeff0->compute();

  shared_ptr<MomentBatch> coeff1;
  if (aux_dec_[num]) {
    coeff1 = shared_ptr<MomentBatch>(new MomentBatch(array<shared_ptr<const Shell>,2>{{shells_[num]->cartesian_shell(), aux_dec_[num]}}, stack_));
    coeff1->compute();
  } else {
    // TODO just to run
    coeff1 = coeff0;
  }

  const double* carea0 = coeff0->data();
  const double* carea1 = coeff1->data();

  shared_ptr<Matrix> tmparea(new Matrix(ssize,a));
  array<shared_ptr<Matrix>,3> out;

  for (int i = 0; i != 3; ++i, carea0 += coeff0->size_block(), carea1 += coeff1->size_block()) {
    if (shells_[num]->spherical()) {
      const int carsphindex = shells_[num]->angular_number() * ANG_HRR_END + 0; // only transform shell
      const int nloop = shells_[num]->num_contracted() * asize_inc; 
      carsphlist.carsphfunc_call(carsphindex, nloop, carea0, tmparea->data());
    } else {
      assert(coeff0->size_block() == asize_inc*ssize);
      copy(carea0, carea0+coeff0->size_block(), tmparea->data());
    }
    if (aux_dec_[num]) {
      if (shells_[num]->spherical()) {
        const int carsphindex = shells_[num]->angular_number() * ANG_HRR_END + 0; // only transform shell
        const int nloop = shells_[num]->num_contracted() * asize_dec; 
        carsphlist.carsphfunc_call(carsphindex, nloop, carea1, tmparea->data()+asize_inc*ssize);
      } else {
        assert(coeff1->size_block() == asize_dec*ssize);
        copy(carea1, carea1+coeff1->size_block(), tmparea->data()+asize_inc*ssize);
        }
    }

    out[i] = shared_ptr<Matrix>(new Matrix(*overlap ^ *tmparea));
  }
  return out;
}

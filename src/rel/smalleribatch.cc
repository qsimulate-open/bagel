//
// BAGEL - Parallel electron correlation program.
// Filename: smalleribatch.cc
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
#include <src/rel/smalleribatch.h>

using namespace std;
using namespace bagel;

SmallERIBatch::SmallERIBatch(std::array<std::shared_ptr<const Shell>,4> info, const double dummy)
  : shells_{{info[1],info[2],info[3]}}, stack_(resources__->get()) {

  assert(info[0]->dummy());

  size_block_ = shells_[0]->nbasis()*shells_[1]->nbasis()*shells_[2]->nbasis();
  size_alloc_ = size_block_ * 4;
  data_ = stack_->get(size_alloc_);
  fill(data_, data_+size_alloc_, 0.0);
}


SmallERIBatch::~SmallERIBatch() {
  stack_->release(size_alloc_, data_);
  resources__->release(stack_);
}

void SmallERIBatch::compute() {

  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int s2size = shells_[2]->nbasis();
  const int a1size_inc = shells_[1]->aux_inc()->nbasis();
  const int a2size_inc = shells_[2]->aux_inc()->nbasis();
  const int a1size_dec = shells_[1]->aux_dec() ? shells_[1]->aux_dec()->nbasis() : 0;
  const int a2size_dec = shells_[2]->aux_dec() ? shells_[2]->aux_dec()->nbasis() : 0;
  const int a1 = a1size_inc + a1size_dec;
  const int a2 = a2size_inc + a2size_dec;

  // first compute uncontracted ERI with auxiliary basis (cartesian)

  double* eri = stack_->get(s0size * a1 * a2);

  eri_compute(eri);

  double* ints = stack_->get(s0size*a1*s2size);

  array<double* const,4> data = {{data_, data_+size_block_, data_+size_block_*2, data_+size_block_*3}};
  array<double,9> coeff = {{1.0,1.0,-1.0,-1.0,1.0,1.0,1.0,-1.0,1.0}};
  array<int,9> b = {{0,1,3,1,0,2,3,2,0}};

  for (int i = 0; i != 3; ++i) {
    dgemm_("N", "N", s0size*a1, s2size, a2, 1.0, eri, s0size*a1, shells_[2]->small(i)->data(), a2, 0.0, ints, s0size*a1);
    for (int k = 0; k != 3; ++k) {
      for (int j = 0; j != s2size; ++j) {
        dgemm_("N", "N", s0size, s1size, a1, coeff[k+3*i], ints+j*s0size*a1, s0size, shells_[1]->small(k)->data(), a1, 1.0, data[b[k+3*i]]+j*s0size*s1size, s0size); 
      }
    }
  }

stack_->release(s0size*a1*s2size, ints);
stack_->release(s0size*a1*a2, eri);
#if 0
  std::array<shared_ptr<Matrix>,3> ints;
  for (int i = 0; i != 3; ++i)
    ints[i] = shared_ptr<Matrix>(new Matrix(*shells_[0]->small(i) % *eri));

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
#endif

}


namespace bagel {
  struct Address {
    const size_t ld0;
    const size_t ld1;
    const size_t ld2;
    Address(const size_t l0, const size_t l1, const size_t l2) : ld0(l0), ld1(l1), ld2(l2) {}
    size_t operator()(const size_t& i, const size_t& j, const size_t& k) const { return i+ld0*(j+ld1*k); }
  };
}


void SmallERIBatch::eri_compute(double* eri) const {

  // shells_[0] is aux function, shelles_[1] and [2] are basis 

  const int s0size = shells_[0]->nbasis();
  const int a1size_inc = shells_[1]->aux_inc()->nbasis();
  const int a2size_inc = shells_[2]->aux_inc()->nbasis();
  const int a1size_dec = shells_[1]->aux_dec() ? shells_[1]->aux_dec()->nbasis() : 0;
  const int a2size_dec = shells_[2]->aux_dec() ? shells_[2]->aux_dec()->nbasis() : 0;
  const int a1 = a1size_inc + a1size_dec;
  const int a2 = a2size_inc + a2size_dec;

  const shared_ptr<const Shell> dummy(new Shell(shells_[0]->spherical()));
  Address m(s0size, a1, a2);

  {
    shared_ptr<ERIBatch> eric(new ERIBatch(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_inc(), shells_[2]->aux_inc()}}, 2.0, 0.0, true, stack_));
    eric->compute();
    for (int i = 0; i != a2size_inc; i++)
      copy_n(eric->data() + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,i));
  }
  if (shells_[1]->aux_dec() && shells_[2]->aux_dec()) {
    shared_ptr<ERIBatch> eric(new ERIBatch(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_dec(), shells_[2]->aux_dec()}}, 2.0, 0.0, true, stack_));
    eric->compute();
    for (int i = 0; i != a2size_dec; i++)
      copy_n(eric->data() + i * s0size * a1size_dec, s0size * a1size_dec, eri + m(0,a1size_inc,a2size_inc+i));
  }
  if (shells_[1]->aux_dec()) {
    shared_ptr<ERIBatch> eric(new ERIBatch(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_dec(), shells_[2]->aux_inc()}}, 2.0, 0.0, true, stack_));
    eric->compute();
    for (int i = 0; i != a2size_inc; i++)
      copy_n(eric->data() + i * s0size * a1size_dec, s0size * a1size_dec, eri + m(0,a1size_inc, i));
  }
  if (shells_[2]->aux_dec()) {
    shared_ptr<ERIBatch> eric(new ERIBatch(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_inc(), shells_[2]->aux_dec()}}, 2.0, 0.0, true, stack_));
    eric->compute();
    for (int i = 0; i != a2size_dec; i++)
      copy_n(eric->data() + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,a2size_inc+i));
  }
}

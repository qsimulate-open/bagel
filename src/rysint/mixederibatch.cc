//
// BAGEL - Parallel electron correlation program.
// Filename: mixederibatch.cc
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
#include <src/rysint/carsphlist.h>
#include <src/rysint/mixederibatch.h>
#include <src/rysint/libint.h>

using namespace std;
using namespace bagel;
const static CarSphList carsphlist;

MixedERIBatch::MixedERIBatch(std::array<std::shared_ptr<const Shell>,4> info, const double dummy)
  : shells_{{info[1],info[2],info[3]}}, stack_(resources__->get()) {

  assert(info[0]->dummy());

  size_block_ = shells_[0]->nbasis()*shells_[1]->nbasis()*shells_[2]->nbasis();
  size_alloc_ = size_block_ * 3;
  data_ = stack_->get(size_alloc_);
}


MixedERIBatch::~MixedERIBatch() {
  stack_->release(size_alloc_, data_);
  resources__->release(stack_);
}

double* MixedERIBatch::data(const int i) {
  return data_+i*size_block_;
}

void MixedERIBatch::compute() {

  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int s2size = shells_[2]->nbasis();
  const int a1size_inc = shells_[1]->aux_inc()->nbasis();
  const int a1size_dec = shells_[1]->aux_dec() ? shells_[1]->aux_dec()->nbasis() : 0;
  const int a1 = a1size_inc + a1size_dec;

  // first compute uncontracted ERI with auxiliary basis (cartesian)

  double* eri = stack_->get(s0size * a1 * s2size);

  eri_compute(eri);

  double* ints = stack_->get(s0size*a1*s2size);

  array<double* const,3> data = {{data_, data_+size_block_, data_+size_block_*2}};

  // We are making 3 blocks, XLarge, YLarge, and ZLarge IN THAT ORDER

  for (int k = 0; k != 3; ++k) {
    for (int j = 0; j != s2size; ++j) {
      dgemm_("N", "N", s0size, s1size, a1, 1.0, eri+j*s0size*a1, s0size, shells_[1]->small(k)->data(), a1, 0.0, data[k]+j*s0size*s1size, s0size); 
    }
  }

  stack_->release(s0size*a1*s2size, ints);
  stack_->release(s0size*a1*s2size, eri);
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


void MixedERIBatch::eri_compute(double* eri) const {

  // shells_[0] is aux function, shelles_[1] and [2] are basis 

  const int s0size = shells_[0]->nbasis();
  const int s2size = shells_[2]->nbasis();
  const int a1size_inc = shells_[1]->aux_inc()->nbasis();
  const int a1size_dec = shells_[1]->aux_dec() ? shells_[1]->aux_dec()->nbasis() : 0;
  const int a1 = a1size_inc + a1size_dec;

  const shared_ptr<const Shell> dummy(new Shell(shells_[0]->spherical()));
  Address m(s0size, a1, s2size);

  {
    shared_ptr<const Shell> cart2 = shells_[2]->cartesian_shell();
    const int s2cart = cart2->nbasis();
#ifndef LIBINT_INTERFACE
    shared_ptr<ERIBatch> eric(new ERIBatch(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_inc(), cart2}}, 2.0, 0.0, true, stack_));
#else
    shared_ptr<Libint> eric(new Libint(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_inc(), cart2}}, 2.0));
#endif
    eric->compute();
    // TODO this could be improved
    double* tmp = stack_->get(s0size * a1size_inc * s2cart);
    double* tmp2 = stack_->get(s0size * a1size_inc * s2size);
    mytranspose_(eric->data(), s0size*a1size_inc, s2cart, tmp);
    const int carsphindex = shells_[2]->angular_number() * ANG_HRR_END;
    carsphlist.carsphfunc_call(carsphindex, s0size*a1size_inc, tmp, tmp2);
    mytranspose_(tmp2, s2size, s0size*a1size_inc, tmp);
    
    for (int i = 0; i != s2size; i++)
      copy_n(tmp + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,i));

    stack_->release(s0size * a1size_inc * s2size, tmp2);
    stack_->release(s0size * a1size_inc * s2cart, tmp);
  }
  if (shells_[1]->aux_dec()) {
    shared_ptr<const Shell> cart2 = shells_[2]->cartesian_shell();
    const int s2cart = cart2->nbasis();
#ifndef LIBINT_INTERFACE 
    shared_ptr<ERIBatch> eric(new ERIBatch(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_dec(), cart2}},
                                           2.0, 0.0, true, stack_));
#else
    shared_ptr<Libint> eric(new Libint(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_dec(), cart2}}, 2.0));
#endif
    eric->compute();
    // TODO this could be improved
    double* tmp = stack_->get(s0size * a1size_dec * s2cart);
    double* tmp2 = stack_->get(s0size * a1size_dec * s2size);
    mytranspose_(eric->data(), s0size*a1size_dec, s2cart, tmp);
    const int carsphindex = shells_[2]->angular_number() * ANG_HRR_END;
    carsphlist.carsphfunc_call(carsphindex, s0size*a1size_dec, tmp, tmp2);
    mytranspose_(tmp2, s2size, s0size*a1size_dec, tmp);

    for (int i = 0; i != s2size; i++)
      copy_n(tmp + i * s0size * a1size_dec, s0size * a1size_dec, eri + m(0, a1size_inc, i));

    stack_->release(s0size * a1size_dec * s2size, tmp2);
    stack_->release(s0size * a1size_dec * s2cart, tmp);
  }
}

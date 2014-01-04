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


#include <src/integral/libint/libint.h>
#include <src/integral/rys/smalleribatch.h>

using namespace std;
using namespace bagel;

SmallERIBatch::SmallERIBatch(std::array<std::shared_ptr<const Shell>,4> info, const double dummy)
  : shells_{{info[1],info[2],info[3]}}, stack_(resources__->get()) {

  assert(info[0]->dummy());

  size_block_ = shells_[0]->nbasis()*shells_[1]->nbasis()*shells_[2]->nbasis();
  size_alloc_ = size_block_ * 6;
  data_ = stack_->get(size_alloc_);
}


SmallERIBatch::~SmallERIBatch() {
  stack_->release(size_alloc_, data_);
  resources__->release(stack_);
}


void SmallERIBatch::compute() {

  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a2size_inc = shells_[2]->nbasis_aux_increment();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  const size_t a2size_dec = shells_[2]->nbasis_aux_decrement();
  const size_t a1 = a1size_inc + a1size_dec;
  const size_t a2 = a2size_inc + a2size_dec;

  const size_t s0size = shells_[0]->nbasis();
  const size_t s1size = shells_[1]->nbasis();
  const size_t s2size = shells_[2]->nbasis();

  // first compute uncontracted ERI with auxiliary basis (cartesian)

  double* eri = stack_->get(s0size * a1 * a2);

  eri_compute(eri);

  double* ints = stack_->get(s0size*a1*s2size);

  array<double* const,6> data = {{data_, data_+size_block_, data_+size_block_*2, data_+size_block_*3, data_+size_block_*4, data_+size_block_*5}};

  for (int i = 0; i != 3; ++i) {
    dgemm_("N", "N", s0size*a1, s2size, a2, 1.0, eri, s0size*a1, shells_[2]->small(i)->data(), a2, 0.0, ints, s0size*a1);
    for (int k = 0; k <= i; ++k) {
      for (int j = 0; j != s2size; ++j) {
        dgemm_("N", "N", s0size, s1size, a1, 1.0, ints+j*s0size*a1, s0size, shells_[1]->small(k)->data(), a1, 0.0, data[k*(5-k)/2+i]+j*s0size*s1size, s0size);
      }
    }
  }

  stack_->release(s0size*a1*s2size, ints);
  stack_->release(s0size*a1*a2, eri);
}


void SmallERIBatch::eri_compute(double* eri) const {

  // shells_[0] is aux function, shelles_[1] and [2] are basis

  const size_t s0size = shells_[0]->nbasis();
  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a2size_inc = shells_[2]->nbasis_aux_increment();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  const size_t a2size_dec = shells_[2]->nbasis_aux_decrement();
  const size_t a1 = a1size_inc + a1size_dec;
  const size_t a2 = a2size_inc + a2size_dec;

  auto dummy = make_shared<const Shell>(shells_[0]->spherical());
  auto m = [&s0size, &a1](const size_t& i, const size_t& j, const size_t k){ return i+s0size*(j+a1*k); }; 

  {
#ifndef LIBINT_INTERFACE
    auto eric = make_shared<ERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_increment()}}, 2.0, 0.0, true, stack_);
#else
    auto eric = make_shared<Libint>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_increment()}}, 0.0, stack_);
#endif
    eric->compute();
    for (int i = 0; i != a2size_inc; i++)
      copy_n(eric->data() + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,i));
  }
  if (shells_[1]->aux_decrement() && shells_[2]->aux_decrement()) {
#ifndef LIBINT_INTERFACE
    auto eric = make_shared<ERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_decrement()}}, 2.0, 0.0, true, stack_);
#else
    auto eric = make_shared<Libint>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_decrement()}}, 0.0, stack_);
#endif
    eric->compute();
    for (int i = 0; i != a2size_dec; i++)
      copy_n(eric->data() + i * s0size * a1size_dec, s0size * a1size_dec, eri + m(0,a1size_inc,a2size_inc+i));
  }
  if (shells_[1]->aux_decrement()) {
#ifndef LIBINT_INTERFACE
    auto eric = make_shared<ERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_increment()}}, 2.0, 0.0, true, stack_);
#else
    auto eric = make_shared<Libint>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_increment()}}, 0.0, stack_);
#endif
    eric->compute();
    for (int i = 0; i != a2size_inc; i++)
      copy_n(eric->data() + i * s0size * a1size_dec, s0size * a1size_dec, eri + m(0,a1size_inc, i));
  }
  if (shells_[2]->aux_decrement()) {
#ifndef LIBINT_INTERFACE
    auto eric = make_shared<ERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_decrement()}}, 2.0, 0.0, true, stack_);
#else
    auto eric = make_shared<Libint>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_decrement()}}, 0.0, stack_);
#endif
    eric->compute();
    for (int i = 0; i != a2size_dec; i++)
      copy_n(eric->data() + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,a2size_inc+i));
  }
}

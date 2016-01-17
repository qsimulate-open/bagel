//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexsmalleribatch.cc
// Copyright (C) 2014 Toru Shiozaki
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


#include <src/integral/comprys/complexsmalleribatch.h>

using namespace std;
using namespace bagel;

ComplexSmallERIBatch::ComplexSmallERIBatch(array<shared_ptr<const Shell>,4> info, const double dummy)
  : shells_{{info[1],info[2],info[3]}}, stack_(resources__->get()) {

  assert(info[0]->dummy());

  size_block_ = shells_[0]->nbasis()*shells_[1]->nbasis()*shells_[2]->nbasis();
  size_alloc_ = size_block_ * 6;
  data_ = stack_->get<complex<double>>(size_alloc_);
}


ComplexSmallERIBatch::~ComplexSmallERIBatch() {
  stack_->release(size_alloc_, data_);
  resources__->release(stack_);
}


void ComplexSmallERIBatch::compute() {

  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a2size_inc = shells_[2]->nbasis_aux_increment();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  const size_t a2size_dec = shells_[2]->nbasis_aux_decrement();
  const size_t a1size_same = shells_[1]->nbasis_aux_same();
  const size_t a2size_same = shells_[2]->nbasis_aux_same();
  const size_t a1 = a1size_inc + a1size_dec + a1size_same;
  const size_t a2 = a2size_inc + a2size_dec + a2size_same;

  const size_t s0size = shells_[0]->nbasis();
  const size_t s1size = shells_[1]->nbasis();
  const size_t s2size = shells_[2]->nbasis();

  // first compute uncontracted ERI with auxiliary basis (cartesian)

  complex<double>* eri = stack_->get<complex<double>>(s0size * a1 * a2);

  eri_compute(eri);

  complex<double>* ints = stack_->get<complex<double>>(s0size*a1*s2size);

  array<complex<double>* const,6> data = {{data_, data_+size_block_, data_+size_block_*2, data_+size_block_*3, data_+size_block_*4, data_+size_block_*5}};

  for (int i = 0; i != 3; ++i) {
    zgemm3m_("N", "N", s0size*a1, s2size, a2, 1.0, eri, s0size*a1, shells_[2]->zsmall(i)->data(), a2, 0.0, ints, s0size*a1);
    for (int k = 0; k <= i; ++k) {
      for (int j = 0; j != s2size; ++j) {
        zgemm3m_("N", "N", s0size, s1size, a1, 1.0, ints+j*s0size*a1, s0size, shells_[1]->zsmallc(k)->data(), a1, 0.0, data[k*(5-k)/2+i]+j*s0size*s1size, s0size);
      }
    }
  }

  stack_->release(s0size*a1*s2size, ints);
  stack_->release(s0size*a1*a2, eri);
}


void ComplexSmallERIBatch::eri_compute(complex<double>* eri) const {

  // shells_[0] is aux function, shells_[1] and [2] are basis

  const size_t s0size = shells_[0]->nbasis();
  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a2size_inc = shells_[2]->nbasis_aux_increment();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  const size_t a2size_dec = shells_[2]->nbasis_aux_decrement();
  const size_t a1size_same = shells_[1]->nbasis_aux_same();
  const size_t a2size_same = shells_[2]->nbasis_aux_same();
  const size_t a1 = a1size_inc + a1size_dec + a1size_same;

  auto dummy = make_shared<const Shell>(shells_[0]->spherical());
  auto m = [&s0size, &a1](const size_t& i, const size_t& j, const size_t k){ return i+s0size*(j+a1*k); };

  {
    auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_increment()}}, 2.0, 0.0, true, stack_);
    eric->compute();
    for (int i = 0; i != a2size_inc; ++i)
      copy_n(eric->data() + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,i));
  }
  if (shells_[1]->aux_decrement() && shells_[2]->aux_decrement()) {
    auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_decrement()}}, 2.0, 0.0, true, stack_);
    eric->compute();
    for (int i = 0; i != a2size_dec; ++i)
      copy_n(eric->data() + i * s0size * a1size_dec, s0size * a1size_dec, eri + m(0,a1size_inc,a2size_inc+i));
  }
  if (shells_[1]->aux_decrement()) {
    auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_increment()}}, 2.0, 0.0, true, stack_);
    eric->compute();
    for (int i = 0; i != a2size_inc; ++i)
      copy_n(eric->data() + i * s0size * a1size_dec, s0size * a1size_dec, eri + m(0,a1size_inc, i));
  }
  if (shells_[2]->aux_decrement()) {
    auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_decrement()}}, 2.0, 0.0, true, stack_);
    eric->compute();
    for (int i = 0; i != a2size_dec; ++i)
      copy_n(eric->data() + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,a2size_inc+i));
  }

  // Unchanged angular momentum (common origin only)
  if (shells_[1]->aux_same()) {
    assert(shells_[2]->aux_same());
    const size_t a1size_id = a1size_inc + a1size_dec;
    const size_t a2size_id = a2size_inc + a2size_dec;
    {
      auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_same()}}, 2.0, 0.0, true, stack_);
      eric->compute();
      for (int i = 0; i != a2size_same; ++i)
        copy_n(eric->data() + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,a2size_id+i));
    }
    {
      auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_same(), shells_[2]->aux_increment()}}, 2.0, 0.0, true, stack_);
      eric->compute();
      for (int i = 0; i != a2size_inc; ++i)
        copy_n(eric->data() + i * s0size * a1size_same, s0size * a1size_same, eri + m(0,a1size_id,i));
    }
    {
      auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_same(), shells_[2]->aux_same()}}, 2.0, 0.0, true, stack_);
      eric->compute();
      for (int i = 0; i != a2size_same; ++i)
        copy_n(eric->data() + i * s0size * a1size_same, s0size * a1size_same, eri + m(0,a1size_id,a2size_id+i));
    }
    if (shells_[1]->aux_decrement()) {
      auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_same()}}, 2.0, 0.0, true, stack_);
      eric->compute();
      for (int i = 0; i != a2size_same; ++i)
        copy_n(eric->data() + i * s0size * a1size_dec, s0size * a1size_dec, eri + m(0,a1size_inc,a2size_id+i));
    }
    if (shells_[2]->aux_decrement()) {
      auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_same(), shells_[2]->aux_decrement()}}, 2.0, 0.0, true, stack_);
      eric->compute();
      for (int i = 0; i != a2size_dec; ++i)
        copy_n(eric->data() + i * s0size * a1size_same, s0size * a1size_same, eri + m(0,a1size_id,a2size_inc+i));
    }
  } else assert(!shells_[2]->aux_same());


}

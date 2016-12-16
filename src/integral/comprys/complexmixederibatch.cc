//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexmixederibatch.cc
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


#include <iomanip>
#include <src/integral/carsphlist.h>
#include <src/integral/comprys/complexmixederibatch.h>

using namespace std;
using namespace bagel;
const static CCarSphList carsphlist;

ComplexMixedERIBatch::ComplexMixedERIBatch(array<shared_ptr<const Shell>,4> info, const double dummy)
  : shells_{{info[1],info[2],info[3]}}, stack_(resources__->get()) {

  assert(info[0]->dummy());

  size_block_ = shells_[0]->nbasis()*shells_[1]->nbasis()*shells_[2]->nbasis();
  size_alloc_ = size_block_ * Nblocks();
  data_ = stack_->get<complex<double>>(size_alloc_);
}


ComplexMixedERIBatch::~ComplexMixedERIBatch() {
  stack_->release(size_alloc_, data_);
  resources__->release(stack_);
}


void ComplexMixedERIBatch::compute() {

  const int s0size = shells_[0]->nbasis();
  const int s1size = shells_[1]->nbasis();
  const int s2size = shells_[2]->nbasis();
  const int a1size_inc = shells_[1]->aux_increment()->nbasis();
  const int a1size_dec = shells_[1]->aux_decrement() ? shells_[1]->aux_decrement()->nbasis() : 0;
  const int a1size_same = shells_[1]->aux_same() ? shells_[1]->aux_same()->nbasis() : 0;
  const int a1 = a1size_inc + a1size_dec + a1size_same;

  // first compute uncontracted ERI with auxiliary basis (cartesian)

  complex<double>* eri = stack_->get<complex<double>>(s0size*a1*s2size);

  eri_compute(eri);

  complex<double>* ints = stack_->get<complex<double>>(s0size*a1*s2size);

  array<complex<double>* const,3> data = {{data_, data_+size_block_, data_+size_block_*2}};

  // We are making 3 blocks, XLarge, YLarge, and ZLarge IN THAT ORDER

  for (int k = 0; k != 3; ++k) {
    for (int j = 0; j != s2size; ++j) {
      zgemm3m_("N", "N", s0size, s1size, a1, 1.0, eri+j*s0size*a1, s0size, shells_[1]->zsmallc(k)->data(), a1, 0.0, data[k]+j*s0size*s1size, s0size);
    }
  }

  stack_->release(s0size*a1*s2size, ints);
  stack_->release(s0size*a1*s2size, eri);
}


void ComplexMixedERIBatch::eri_compute(complex<double>* eri) const {

  // shells_[0] is aux function, shells_[1] and [2] are basis

  const int s0size = shells_[0]->nbasis();
  const int s2size = shells_[2]->nbasis();
  const int a1size_inc = shells_[1]->aux_increment()->nbasis();
  const int a1size_dec = shells_[1]->aux_decrement() ? shells_[1]->aux_decrement()->nbasis() : 0;
  const int a1size_same = shells_[1]->aux_same() ? shells_[1]->aux_same()->nbasis() : 0;
  const int a1 = a1size_inc + a1size_dec + a1size_same;

  auto dummy = make_shared<const Shell>(shells_[0]->spherical());
  auto m = [&s0size, &a1](const size_t& i, const size_t& j, const size_t k){ return i+s0size*(j+a1*k); };

  {
    shared_ptr<const Shell> cart2 = shells_[2]->cartesian_shell();
    const int s2cart = cart2->nbasis();
    auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), cart2}}, 2.0, 0.0, true, stack_);
    eric->compute();

    complex<double>* tmp = stack_->get<complex<double>>(s0size*a1size_inc*s2cart);
    if (shells_[1]->spherical()) {
      // TODO this could be improved
      complex<double>* tmp2 = stack_->get<complex<double>>(s0size*a1size_inc*s2size);
      blas::transpose(eric->data(), s0size*a1size_inc, s2cart, tmp);

      const int carsphindex = shells_[2]->angular_number() * ANG_HRR_END;
      carsphlist.carsphfunc_call(carsphindex, s0size*a1size_inc*cart2->num_contracted(), tmp, tmp2);

      blas::transpose(tmp2, s2size, s0size*a1size_inc, tmp);
      stack_->release(s0size*a1size_inc*s2size, tmp2);
    } else {
      copy_n(eric->data(), s0size*a1size_inc*s2cart, tmp);
    }

    for (int i = 0; i != s2size; ++i)
      copy_n(tmp + i * s0size * a1size_inc, s0size * a1size_inc, eri + m(0,0,i));

    stack_->release(s0size * a1size_inc * s2cart, tmp);
  }
  if (shells_[1]->aux_decrement()) {
    shared_ptr<const Shell> cart2 = shells_[2]->cartesian_shell();
    const int s2cart = cart2->nbasis();
    auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), cart2}},
                                      2.0, 0.0, true, stack_);
    eric->compute();

    complex<double>* tmp = stack_->get<complex<double>>(s0size*a1size_dec*s2cart);
    if (shells_[1]->spherical()) {
      // TODO this could be improved
      complex<double>* tmp2 = stack_->get<complex<double>>(s0size*a1size_dec*s2size);
      blas::transpose(eric->data(), s0size*a1size_dec, s2cart, tmp);

      const int carsphindex = shells_[2]->angular_number() * ANG_HRR_END;
      carsphlist.carsphfunc_call(carsphindex, s0size*a1size_dec*cart2->num_contracted(), tmp, tmp2);

      blas::transpose(tmp2, s2size, s0size*a1size_dec, tmp);
      stack_->release(s0size*a1size_dec*s2size, tmp2);
    } else {
      copy_n(eric->data(), s0size*a1size_dec*s2cart, tmp);
    }

    for (int i = 0; i != s2size; ++i)
      copy_n(tmp + i*s0size*a1size_dec, s0size*a1size_dec, eri + m(0,a1size_inc,i));

    stack_->release(s0size * a1size_dec * s2cart, tmp);
  }

  // Unchanged angular momentum (common origin only)
  if (shells_[1]->aux_same()) {
    const size_t a1size_id = a1size_inc + a1size_dec;
    shared_ptr<const Shell> cart2 = shells_[2]->cartesian_shell();
    const int s2cart = cart2->nbasis();
    auto eric = make_shared<ComplexERIBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_same(), cart2}},
                                      2.0, 0.0, true, stack_);
    eric->compute();

    complex<double>* tmp = stack_->get<complex<double>>(s0size*a1size_same*s2cart);
    if (shells_[1]->spherical()) {
      // TODO this could be improved
      complex<double>* tmp2 = stack_->get<complex<double>>(s0size*a1size_same*s2size);
      blas::transpose(eric->data(), s0size*a1size_same, s2cart, tmp);

      const int carsphindex = shells_[2]->angular_number() * ANG_HRR_END;
      carsphlist.carsphfunc_call(carsphindex, s0size*a1size_same*cart2->num_contracted(), tmp, tmp2);

      blas::transpose(tmp2, s2size, s0size*a1size_same, tmp);
      stack_->release(s0size*a1size_same*s2size, tmp2);
    } else {
      copy_n(eric->data(), s0size*a1size_same*s2cart, tmp);
    }

    for (int i = 0; i != s2size; ++i)
      copy_n(tmp + i*s0size*a1size_same, s0size*a1size_same, eri + m(0,a1size_id,i));

    stack_->release(s0size * a1size_same * s2cart, tmp);
  }



}

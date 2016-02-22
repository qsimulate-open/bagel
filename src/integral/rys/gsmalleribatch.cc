//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gsmalleribatch.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/util/alpha.h>
#include <src/integral/rys/gradbatch.h>
#include <src/integral/libint/glibint.h>
#include <src/integral/rys/gsmalleribatch.h>

using namespace std;
using namespace bagel;

GSmallERIBatch::GSmallERIBatch(array<shared_ptr<const Shell>,4> info, array<int,3> atoms, const int natoms)
 : shells_{{info[1],info[2],info[3]}}, atoms_(atoms), natoms_(natoms), stack_(resources__->get()) {

  // shells_[0] is aux function, shelles_[1] and [2] are basis
  const size_t s0size = shells_[0]->nbasis();
  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a2size_inc = shells_[2]->nbasis_aux_increment();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  const size_t a2size_dec = shells_[2]->nbasis_aux_decrement();

  assert(info[0]->dummy());

  // we store primitive integrals
  size_block_ = s0size*(a1size_inc + a1size_dec)*(a2size_inc + a2size_dec);
  size_alloc_ = size_block_*9;
  data_ = stack_->get(size_alloc_);
}


GSmallERIBatch::~GSmallERIBatch() {
  stack_->release(size_alloc_, data_);
  resources__->release(stack_);
}


shared_ptr<GradFile> GSmallERIBatch::compute_gradient(array<shared_ptr<const btas::Tensor3<double>>,6>& d) const {
  auto out = make_shared<GradFile>(natoms_);
  static_assert(Comp::X == 0 && Comp::Y == 1 && Comp::Z == 2, "something is wrong in GSmallERIBatch::compute_gradient");

  const size_t a1 = shells_[1]->nbasis_aux_increment() + shells_[1]->nbasis_aux_decrement();
  const size_t a2 = shells_[2]->nbasis_aux_increment() + shells_[2]->nbasis_aux_decrement();

  const size_t s0size = shells_[0]->nbasis();
  const size_t s1size = shells_[1]->nbasis();
  const size_t s2size = shells_[2]->nbasis();

  double* denc = stack_->get(s0size*a1*a2);
  fill_n(denc, s0size*a1*a2, 0.0);

  double* tmp = stack_->get(s0size*s1size*a2);

  int cnt = 0;
  array<int,3> xyz{{Comp::X, Comp::Y, Comp::Z}};
  for (int& i : xyz) {
    for (int& j : xyz) {
      if (i <= j) {
        dgemm_("N", "T", s0size*s1size, a2, s2size, 1.0, d[cnt++]->data(), s0size*s1size, shells_[2]->small(j)->data(), a2, 0.0, tmp, s0size*s1size);
        for (int a = 0; a != a2; ++a) {
          dgemm_("N", "T", s0size, a1, s1size, 1.0, tmp+a*s0size*s1size, s0size, shells_[1]->small(i)->data(), a1, 1.0, denc+a*s0size*a1, s0size);
        }
      }
    }
  }

  for (int iatom = 0; iatom != 3; ++iatom) {
    for (int i = 0; i != 3; ++i) {
      out->element(i, atoms_[iatom]) += ddot_(s0size*a1*a2, denc, 1, data(i+3*iatom), 1);
    }
  }

  stack_->release(s0size*s1size*a2, tmp);
  stack_->release(s0size*a1*a2, denc);
  return out;
}


void GSmallERIBatch::compute() {

  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a2size_inc = shells_[2]->nbasis_aux_increment();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  const size_t a2size_dec = shells_[2]->nbasis_aux_decrement();
  const size_t a1 = a1size_inc + a1size_dec;

  auto dummy = make_shared<const Shell>(shells_[0]->spherical());
  const size_t s0size = shells_[0]->nbasis();
  auto m = [&s0size, &a1](const size_t& i, const size_t& j, const size_t k){ return i+s0size*(j+a1*k); };

  {
#ifndef LIBINT_INTERFACE
    auto eric = make_shared<GradBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_increment()}}, 2.0, 0.0, true, stack_);
#else
    auto eric = make_shared<GLibint>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_increment()}}, stack_);
#endif
    eric->compute();

    array<int,4> jatom = {{-1, 2, 1, 0}};
    if (eric->swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
    if (eric->swap01()) swap(jatom[0], jatom[1]);
    if (eric->swap23()) swap(jatom[2], jatom[3]);
    for (int iatom = 0; iatom != 4; ++iatom) {
      if (jatom[iatom] < 0) continue;
      for (int icart = 0; icart != 3; ++icart) {
        const int sblock = icart + iatom*3;
        const int iblock = icart + jatom[iatom]*3;
        for (int i = 0; i != a2size_inc; i++)
          copy_n(eric->data(sblock) + i*s0size*a1size_inc, s0size*a1size_inc, data(iblock) + m(0,0,i));
      }
    }
  }
  if (shells_[1]->aux_decrement() && shells_[2]->aux_decrement()) {
#ifndef LIBINT_INTERFACE
    auto eric = make_shared<GradBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_decrement()}}, 2.0, 0.0, true, stack_);
#else
    auto eric = make_shared<GLibint>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_decrement()}}, stack_);
#endif
    eric->compute();

    array<int,4> jatom = {{-1, 2, 1, 0}};
    if (eric->swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
    if (eric->swap01()) swap(jatom[0], jatom[1]);
    if (eric->swap23()) swap(jatom[2], jatom[3]);
    for (int iatom = 0; iatom != 4; ++iatom) {
      if (jatom[iatom] < 0) continue;
      for (int icart = 0; icart != 3; ++icart) {
        const int sblock = icart + iatom*3;
        const int iblock = icart + jatom[iatom]*3;
        for (int i = 0; i != a2size_dec; i++)
          copy_n(eric->data(sblock) + i*s0size*a1size_dec, s0size*a1size_dec, data(iblock) + m(0,a1size_inc,a2size_inc+i));
      }
    }
  }
  if (shells_[1]->aux_decrement()) {
#ifndef LIBINT_INTERFACE
    auto eric = make_shared<GradBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_increment()}}, 2.0, 0.0, true, stack_);
#else
    auto eric = make_shared<GLibint>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_decrement(), shells_[2]->aux_increment()}}, stack_);
#endif
    eric->compute();

    array<int,4> jatom = {{-1, 2, 1, 0}};
    if (eric->swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
    if (eric->swap01()) swap(jatom[0], jatom[1]);
    if (eric->swap23()) swap(jatom[2], jatom[3]);
    for (int iatom = 0; iatom != 4; ++iatom) {
      if (jatom[iatom] < 0) continue;
      for (int icart = 0; icart != 3; ++icart) {
        const int sblock = icart + iatom*3;
        const int iblock = icart + jatom[iatom]*3;
        for (int i = 0; i != a2size_inc; i++)
          copy_n(eric->data(sblock) + i*s0size*a1size_dec, s0size*a1size_dec, data(iblock) + m(0,a1size_inc, i));
      }
    }
  }
  if (shells_[2]->aux_decrement()) {
#ifndef LIBINT_INTERFACE
    auto eric = make_shared<GradBatch>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_decrement()}}, 2.0, 0.0, true, stack_);
#else
    auto eric = make_shared<GLibint>(array<shared_ptr<const Shell>,4>{{dummy, shells_[0], shells_[1]->aux_increment(), shells_[2]->aux_decrement()}}, stack_);
#endif
    eric->compute();

    array<int,4> jatom = {{-1, 2, 1, 0}};
    if (eric->swap0123()) { swap(jatom[0], jatom[2]); swap(jatom[1], jatom[3]); }
    if (eric->swap01()) swap(jatom[0], jatom[1]);
    if (eric->swap23()) swap(jatom[2], jatom[3]);
    for (int iatom = 0; iatom != 4; ++iatom) {
      if (jatom[iatom] < 0) continue;
      for (int icart = 0; icart != 3; ++icart) {
        const int sblock = icart + iatom*3;
        const int iblock = icart + jatom[iatom]*3;
        for (int i = 0; i != a2size_dec; i++)
          copy_n(eric->data(sblock) + i*s0size*a1size_inc, s0size*a1size_inc, data(iblock) + m(0,0,a2size_inc+i));
      }
    }
  }
}

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gsmallnaibatch.cc
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
#include <src/integral/rys/gnaibatch.h>
#include <src/integral/rys/gsmallnaibatch.h>

using namespace std;
using namespace bagel;

GSmallNAIBatch::GSmallNAIBatch(array<shared_ptr<const Shell>,2> info, shared_ptr<const Molecule> mol, const tuple<int,int> i)
    : mol_(mol), shells_(info), size_block_(shells_[0]->nbasis() * shells_[1]->nbasis()), iatom_(i) {

  assert(shells_[0]->relativistic() && shells_[1]->relativistic());
  const size_t a0size_inc = shells_[0]->nbasis_aux_increment();
  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a0size_dec = shells_[0]->nbasis_aux_decrement();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  const size_t a0 = a0size_inc + a0size_dec;
  const size_t a1 = a1size_inc + a1size_dec;

  for (int i = 0; i != mol_->natom()*3; ++i)
    data_.push_back(make_shared<Matrix>(a0, a1, true));
}


shared_ptr<GradFile> GSmallNAIBatch::compute_gradient(array<shared_ptr<const Matrix>,6> d) const {
  auto out = make_shared<GradFile>(mol_->natom());
  static_assert(Comp::X == 0 && Comp::Y == 1 && Comp::Z == 2, "something is wrong in GSmallNAIBatch::compute_gradient");
  array<int,3> xyz{{Comp::X, Comp::Y, Comp::Z}};

  const size_t a0size_inc = shells_[0]->nbasis_aux_increment();
  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a0size_dec = shells_[0]->nbasis_aux_decrement();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  const size_t a0 = a0size_inc + a0size_dec;
  const size_t a1 = a1size_inc + a1size_dec;

  Matrix denc(a0, a1, true);
  int cnt = 0;
  for (int& i : xyz)
    for (int& j : xyz)
      if (i <= j) {
        if (d[cnt])
          denc += *shells_[0]->small(i) * *d[cnt] ^ *shells_[1]->small(j);
        ++cnt;
      }

  for (int iatom = 0; iatom != mol_->natom(); ++iatom)
    for (int i = 0; i != 3; ++i)
      out->element(i, iatom) += data_[i+3*iatom]->dot_product(denc);
  return out;
}


void GSmallNAIBatch::compute() {
  const size_t a0size_inc = shells_[0]->nbasis_aux_increment();
  const size_t a1size_inc = shells_[1]->nbasis_aux_increment();
  const size_t a0size_dec = shells_[0]->nbasis_aux_decrement();
  const size_t a1size_dec = shells_[1]->nbasis_aux_decrement();
  {
    auto naic = make_shared<GNAIBatch>(array<shared_ptr<const Shell>,2>{{shells_[0]->aux_increment(), shells_[1]->aux_increment()}}, mol_, iatom_);
    naic->compute();
    assert(naic->nblocks() == mol_->natom()*3);
    for (int i = 0; i != naic->nblocks(); ++i)
      data_[i]->copy_block(0, 0, a0size_inc, a1size_inc, naic->data(i));
  }
  if (shells_[0]->aux_decrement() && shells_[1]->aux_decrement()) {
    auto naic = make_shared<GNAIBatch>(array<shared_ptr<const Shell>,2>{{shells_[0]->aux_decrement(), shells_[1]->aux_decrement()}}, mol_, iatom_);
    naic->compute();
    for (int i = 0; i != naic->nblocks(); ++i)
      data_[i]->copy_block(a0size_inc, a1size_inc, a0size_dec, a1size_dec, naic->data(i));
  }
  if (shells_[0]->aux_decrement()) {
    auto naic = make_shared<GNAIBatch>(array<shared_ptr<const Shell>,2>{{shells_[0]->aux_decrement(), shells_[1]->aux_increment()}}, mol_, iatom_);
    naic->compute();
    for (int i = 0; i != naic->nblocks(); ++i)
      data_[i]->copy_block(a0size_inc, 0, a0size_dec, a1size_inc, naic->data(i));
  }
  if (shells_[1]->aux_decrement()) {
    auto naic = make_shared<GNAIBatch>(array<shared_ptr<const Shell>,2>{{shells_[0]->aux_increment(), shells_[1]->aux_decrement()}}, mol_, iatom_);
    naic->compute();
    for (int i = 0; i != naic->nblocks(); ++i)
      data_[i]->copy_block(0, a1size_inc, a0size_inc, a1size_dec, naic->data(i));
  }
}

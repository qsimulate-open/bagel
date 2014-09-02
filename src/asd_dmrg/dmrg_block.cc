//
// BAGEL - Parallel electron correlation program.
// Filename: dmrg_block.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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

#include <algorithm>

#include <src/asd_dmrg/dmrg_block.h>
#include <src/ras/civector.h>

using namespace bagel;
using namespace std;

DMRG_Block1::DMRG_Block1(GammaForestASD<RASDvec>&& forest, const map<BlockKey, shared_ptr<const Matrix>> h2e, const map<BlockKey,
           shared_ptr<const Matrix>> spin, shared_ptr<const Matrix> coeff) : DMRG_Block(coeff), H2e_(h2e), spin_(spin) {
  // Build set of blocks
  for (auto& i : h2e) {
    assert(i.second->ndim() == i.second->mdim());
    blocks_.emplace(i.first.nelea, i.first.neleb, i.second->ndim());
  }

  // initialize operator space
  for (auto o : forest.sparselist()) {
    list<GammaSQ> gammalist = get<0>(o);

    BlockInfo brakey = get<1>(o);
    BlockInfo ketkey = get<2>(o);

    assert(blocks_.count(brakey));
    assert(blocks_.count(ketkey));

    const int bratag = forest.block_tag(brakey);
    const int kettag = forest.block_tag(ketkey);

    assert(forest.template exist<0>(bratag, kettag, gammalist));
    shared_ptr<const Matrix> mat = forest.template get<0>(bratag, kettag, gammalist);
    btas::CRange<3> range(brakey.nstates, ketkey.nstates, lrint(pow(norb(), gammalist.size())));
    // checking ndim
    assert(mat->ndim() == range.extent(0)*range.extent(1));
    // checking mdim
    assert(mat->mdim() == range.extent(2));
    // internal check
    assert(mat->storage().size() == range.area());
    // convert this matrix to 3-tensor
    auto tensor = make_shared<btas::Tensor3<double>>(range, move(mat->storage()));
    // add matrix
    CouplingBlock cb(brakey, ketkey, tensor);
    sparse_[gammalist].emplace(cb.key(), cb);
  }
}

DMRG_Block1::DMRG_Block1(GammaForestProdASD&& forest, const map<BlockKey, shared_ptr<const Matrix>> h2e,
                                                    const map<BlockKey, shared_ptr<const Matrix>> spin,
                                                    shared_ptr<const Matrix> coeff) : DMRG_Block(coeff), H2e_(h2e), spin_(spin) {
  // Build set of blocks
  for (auto& i : h2e) {
    assert(i.second->ndim() == i.second->mdim());
    blocks_.emplace(i.first.nelea, i.first.neleb, i.second->ndim());
  }

  // initialize operator space
  for (auto o : forest.sparselist()) {
    list<GammaSQ> gammalist = get<0>(o);

    BlockInfo brakey = get<1>(o);
    BlockInfo ketkey = get<2>(o);

    assert(blocks_.count(brakey));
    assert(blocks_.count(ketkey));

    assert(forest.exist(brakey, ketkey, gammalist));
    shared_ptr<const Matrix> mat = forest.get(brakey, ketkey, gammalist);
    btas::CRange<3> range(brakey.nstates, ketkey.nstates, lrint(pow(norb(), gammalist.size())));
    // checking ndim
    assert(mat->ndim() == range.extent(0)*range.extent(1));
    // checking mdim
    assert(mat->mdim() == range.extent(2));
    // internal check
    assert(mat->storage().size() == range.area());
    // convert this matrix to 3-tensor
    auto tensor = make_shared<btas::Tensor3<double>>(range, move(mat->storage()));
    // add matrix
    CouplingBlock cb(brakey, ketkey, tensor);
    sparse_[gammalist].emplace(cb.key(), cb);
  }
}

shared_ptr<Matrix> DMRG_Block1::spin_lower(const BlockKey k) const {
  BlockKey lowered_key(k.nelea-1, k.neleb+1);
  assert(contains(lowered_key));
  shared_ptr<const btas::Tensor3<double>> gamma = coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}).at({lowered_key, k}).data;

  auto out = make_shared<Matrix>(gamma->extent(0), gamma->extent(1));

  for (int r = 0; r < norb(); ++r)
    blas::ax_plus_y_n(1.0, &(*gamma)(0, 0, r + r*norb()), out->size(), out->data());

  return out;
}

shared_ptr<Matrix> DMRG_Block1::spin_raise(const BlockKey k) const {
  BlockKey raised_key(k.nelea+1, k.neleb-1);
  assert(contains(raised_key));
  shared_ptr<const btas::Tensor3<double>> gamma = coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}).at({k, raised_key}).data;

  auto out = make_shared<Matrix>(gamma->extent(0), gamma->extent(1));

  for (int r = 0; r < norb(); ++r)
    blas::ax_plus_y_n(1.0, &(*gamma)(0, 0, r + r*norb()), out->size(), out->data());

  return out->transpose();
}

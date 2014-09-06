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
#include <src/asd_dmrg/kronecker.h>

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

shared_ptr<const BlockOperators> DMRG_Block1::compute_block_ops(shared_ptr<DimerJop> jop) const {
  return make_shared<BlockOperators1>(shared_from_this(), jop);
}

DMRG_Block2::DMRG_Block2(shared_ptr<const DMRG_Block1> lb, std::shared_ptr<const DMRG_Block1> rb) : left_block_(lb), right_block_(rb) {
  // left block runs first in the resulting pairmap_ vectors
  for (auto& left : lb->blocks()) {
    for (auto& right : rb->blocks()) {
      BlockKey bk(left.nelea+right.nelea, left.neleb+right.neleb);
      const int offset = pairmap_.find(bk)==pairmap_.end() ? 0 : pairmap_[bk].back().offset + pairmap_[bk].back().nstates();
      pairmap_[bk].emplace_back(left, right, offset);
    }
  }

  for (auto& p : pairmap_) {
    BlockKey key = p.first;
    const int nstates = accumulate(p.second.begin(), p.second.end(), 0, [] (int x, DMRG::BlockPair bp) { return x + bp.nstates(); });
    blocks_.insert(BlockInfo(key.nelea, key.neleb, nstates));
  }

  coeff_ = lb->coeff()->merge(rb->coeff());
}

shared_ptr<Matrix> DMRG_Block2::spin(const BlockKey b) const {
  const vector<DMRG::BlockPair>& sec_pairs = pairmap_.at(b);
  BlockInfo binfo = blockinfo(b);
  auto out = make_shared<Matrix>(binfo.nstates, binfo.nstates);

  for (auto& bp : sec_pairs) {
    { // S^2_L
      auto lspin = left_block_->spin(bp.left);
      Matrix eye(bp.right.nstates, bp.right.nstates);
      eye.unit();
      out->add_block(1.0, bp.offset, bp.offset, bp.nstates(), bp.nstates(), kronecker_product(false, eye, false, *lspin));
    }

    { // S^2_R
      auto rspin = right_block_->spin(bp.right);
      Matrix eye(bp.left.nstates, bp.left.nstates);
      eye.unit();
      out->add_block(1.0, bp.offset, bp.offset, bp.nstates(), bp.nstates(), kronecker_product(false, *rspin, false, eye));
    }

    { // 2 * S^z_A S^z_B
      const double sza_szb = 0.5 * static_cast<double>((bp.left.nelea-bp.left.neleb)*(bp.right.nelea-bp.right.neleb));
      for (int i = 0; i < bp.nstates(); ++i)
        out->element(i+bp.offset, i+bp.offset) += sza_szb;
    }

    { // S^+_L S^-_R
      BlockKey lk(bp.left.nelea+1, bp.left.neleb-1);
      BlockKey rk(bp.right.nelea-1, bp.right.neleb+1);

      auto iter = find_if(sec_pairs.begin(), sec_pairs.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk, rk); }
      );

      if (iter != sec_pairs.end()) {
        DMRG::BlockPair tp = *iter;
        auto lowerL = left_block_->spin_lower(bp.left);
        auto raiseR = right_block_->spin_raise(bp.right);

        out->add_block(1.0, tp.offset, bp.offset, tp.nstates(), bp.nstates(), kronecker_product(false, *raiseR, false, *lowerL));
      }
    }

    { // S^-_L S^+_R
      BlockKey lk(bp.left.nelea-1, bp.left.neleb+1);
      BlockKey rk(bp.right.nelea+1, bp.right.neleb-1);

      auto iter = find_if(sec_pairs.begin(), sec_pairs.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk, rk); }
      );

      if (iter != sec_pairs.end()) {
        DMRG::BlockPair tp = *iter;
        auto raiseL = left_block_->spin_raise(bp.left);
        auto lowerR = right_block_->spin_lower(bp.right);

        out->add_block(1.0, tp.offset, bp.offset, tp.nstates(), bp.nstates(), kronecker_product(false, *lowerR, false, *raiseL));
      }
    }
  }

  return out;
}

shared_ptr<Matrix> DMRG_Block2::spin_lower(const BlockKey b) const {
  assert(contains(BlockKey(b.nelea-1, b.neleb+1)));

  BlockInfo source_info = blockinfo(b);
  BlockInfo target_info = blockinfo(BlockKey(b.nelea-1, b.neleb+1));

  auto out = make_shared<Matrix>(target_info.nstates, source_info.nstates);

  const vector<DMRG::BlockPair> pvec = pairmap_.at(b);
  for (auto& source : pvec) {
    { // lower on the left block
      BlockKey lk(source.left.nelea-1, source.left.neleb+1);
      BlockKey rk = source.right.key();
      auto iter = find_if(pvec.begin(), pvec.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk,rk);
      });
      if (iter != pvec.end()) {
        DMRG::BlockPair tp = *iter;
        auto lowerL = left_block_->spin_lower(source.left);
        Matrix eyeR(source.right.nstates, source.right.nstates); eyeR.unit();

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, eyeR, false, *lowerL));
      }
    }

    { // lower on the right block
      BlockKey lk = source.left.key();
      BlockKey rk(source.right.nelea-1, source.right.neleb+1);
      auto iter = find_if(pvec.begin(), pvec.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk,rk);
      });
      if (iter != pvec.end()) {
        DMRG::BlockPair tp = *iter;
        Matrix eyeL(source.left.nstates, source.left.nstates); eyeL.unit();
        auto lowerR = right_block_->spin_lower(source.right);

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, *lowerR, false, eyeL));
      }
    }
  }

  return out;
}

shared_ptr<Matrix> DMRG_Block2::spin_raise(const BlockKey b) const {
  assert(contains(BlockKey(b.nelea+1, b.neleb-1)));

  BlockInfo source_info = blockinfo(b);
  BlockInfo target_info = blockinfo(BlockKey(b.nelea+1, b.neleb-1));

  auto out = make_shared<Matrix>(target_info.nstates, source_info.nstates);

  const vector<DMRG::BlockPair> pvec = pairmap_.at(b);
  for (auto& source : pvec) {
    { // raise on the left block
      BlockKey lk(source.left.nelea+1, source.left.neleb-1);
      BlockKey rk = source.right.key();
      auto iter = find_if(pvec.begin(), pvec.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk,rk);
      });
      if (iter != pvec.end()) {
        DMRG::BlockPair tp = *iter;
        auto raiseL = left_block_->spin_raise(source.left);
        Matrix eyeR(source.right.nstates, source.right.nstates); eyeR.unit();

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, eyeR, false, *raiseL));
      }
    }

    { // raise on the right block
      BlockKey lk = source.left.key();
      BlockKey rk(source.right.nelea+1, source.right.neleb-1);
      auto iter = find_if(pvec.begin(), pvec.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk,rk);
      });
      if (iter != pvec.end()) {
        DMRG::BlockPair tp = *iter;
        Matrix eyeL(source.left.nstates, source.left.nstates); eyeL.unit();
        auto raiseR = right_block_->spin_raise(source.right);

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, *raiseR, false, eyeL));
      }
    }
  }

  return out;
}

shared_ptr<const BlockOperators> DMRG_Block2::compute_block_ops(std::shared_ptr<DimerJop> jop) const {
  return make_shared<BlockOperators2>(shared_from_this(), jop);
}

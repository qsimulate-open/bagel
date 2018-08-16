//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dmrg_block.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <algorithm>

#include <src/asd/dmrg/dmrg_block.h>
#include <src/asd/dmrg/kronecker.h>

using namespace bagel;
using namespace std;

string DMRG_Block1::block_info_to_string(const BlockKey bk, const int state) const {
  stringstream ss;
  ss << "|na:" << bk.nelea << ",nb:" << bk.neleb << "," << state << ">";
  return ss.str();
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

DMRG_Block2::DMRG_Block2(shared_ptr<const DMRG_Block1> lb, shared_ptr<const DMRG_Block1> rb) : left_block_(lb), right_block_(rb) {
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

string DMRG_Block2::block_info_to_string(const BlockKey bk, const int state) const {
  stringstream ss;
  for (auto& bp : pairmap_.at(bk)) {
    if (state >= bp.offset && state < bp.offset+bp.nstates()) {
      const int lstate = (state - bp.offset) % bp.left.nstates;
      const int rstate = (state - bp.offset) / bp.left.nstates;
      ss << left_block_->block_info_to_string(bp.left, lstate) << " (x) " << right_block_->block_info_to_string(bp.right, rstate);
      break;
    }
  }
  return ss.str();
}

shared_ptr<Matrix> DMRG_Block2::spin(const BlockKey b) const {
  const vector<DMRG::BlockPair>& sec_pairs = pairmap_.at(b);
  BlockInfo binfo = blockinfo(b);
  auto out = make_shared<Matrix>(binfo.nstates, binfo.nstates);

  for (auto& source : sec_pairs) {
    { // S^2_L
      auto lspin = left_block_->spin(source.left);
      Matrix eye(source.right.nstates, source.right.nstates);
      eye.unit();
      out->add_block(1.0, source.offset, source.offset, source.nstates(), source.nstates(), kronecker_product(false, eye, false, *lspin));
    }

    { // S^2_R
      auto rspin = right_block_->spin(source.right);
      Matrix eye(source.left.nstates, source.left.nstates);
      eye.unit();
      out->add_block(1.0, source.offset, source.offset, source.nstates(), source.nstates(), kronecker_product(false, *rspin, false, eye));
    }

    { // 2 * S^z_A S^z_B
      auto rspin = right_block_->spin(source.right);
      const double sza_szb = 0.5 * static_cast<double>((source.left.nelea-source.left.neleb)*(source.right.nelea-source.right.neleb));
      for (int i = 0; i < source.nstates(); ++i)
        out->element(i+source.offset, i+source.offset) += sza_szb;
    }

    { // S^-_L S^+_R
      const BlockKey lk(source.left.nelea-1, source.left.neleb+1);
      const BlockKey rk(source.right.nelea+1, source.right.neleb-1);

      auto iter = find_if(sec_pairs.begin(), sec_pairs.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk, rk); }
      );

      if (iter != sec_pairs.end()) {
        const DMRG::BlockPair tp = *iter;
        auto lowerL = left_block_->spin_lower(source.left);
        auto raiseR = right_block_->spin_raise(source.right);

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, *raiseR, false, *lowerL));
      }
    }

    { // S^+_L S^-_R
      const BlockKey lk(source.left.nelea+1, source.left.neleb-1);
      const BlockKey rk(source.right.nelea-1, source.right.neleb+1);

      auto iter = find_if(sec_pairs.begin(), sec_pairs.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk, rk); }
      );

      if (iter != sec_pairs.end()) {
        DMRG::BlockPair tp = *iter;
        auto raiseL = left_block_->spin_raise(source.left);
        auto lowerR = right_block_->spin_lower(source.right);

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, *lowerR, false, *raiseL));
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

  const vector<DMRG::BlockPair> target_pvec = pairmap_.at(target_info.key());
  for (auto& source : pairmap_.at(b)) {
    { // lower on the left block
      BlockKey lk(source.left.nelea-1, source.left.neleb+1);
      BlockKey rk = source.right.key();
      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk,rk);
      });
      if (iter != target_pvec.end()) {
        DMRG::BlockPair tp = *iter;
        auto lowerL = left_block_->spin_lower(source.left);
        Matrix eyeR(source.right.nstates, source.right.nstates); eyeR.unit();

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, eyeR, false, *lowerL));
      }
    }

    { // lower on the right block
      BlockKey lk = source.left.key();
      BlockKey rk(source.right.nelea-1, source.right.neleb+1);
      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk,rk);
      });
      if (iter != target_pvec.end()) {
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

  const vector<DMRG::BlockPair> target_pvec = pairmap_.at(target_info.key());
  for (auto& source : pairmap_.at(b)) {
    { // raise on the left block
      BlockKey lk(source.left.nelea+1, source.left.neleb-1);
      BlockKey rk = source.right.key();
      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk,rk);
      });
      if (iter != target_pvec.end()) {
        DMRG::BlockPair tp = *iter;
        auto raiseL = left_block_->spin_raise(source.left);
        Matrix eyeR(source.right.nstates, source.right.nstates); eyeR.unit();

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, eyeR, false, *raiseL));
      }
    }

    { // raise on the right block
      BlockKey lk = source.left.key();
      BlockKey rk(source.right.nelea+1, source.right.neleb-1);
      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&lk, &rk] (const DMRG::BlockPair& bp) {
        return make_pair(bp.left.key(), bp.right.key())==make_pair(lk,rk);
      });
      if (iter != target_pvec.end()) {
        DMRG::BlockPair tp = *iter;
        Matrix eyeL(source.left.nstates, source.left.nstates); eyeL.unit();
        auto raiseR = right_block_->spin_raise(source.right);

        out->add_block(1.0, tp.offset, source.offset, tp.nstates(), source.nstates(), kronecker_product(false, *raiseR, false, eyeL));
      }
    }
  }

  return out;
}

shared_ptr<const BlockOperators> DMRG_Block2::compute_block_ops(shared_ptr<DimerJop> jop) const {
  return make_shared<BlockOperators2>(shared_from_this(), jop);
}

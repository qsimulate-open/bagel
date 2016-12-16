//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg/block_ops_2.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/dmrg/block_operators.h>
#include <src/asd/dmrg/dmrg_block.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;

BlockOperators2::BlockOperators2(shared_ptr<const DMRG_Block2> blocks, shared_ptr<DimerJop> jop, const double thresh) : BlockOperators(thresh), blocks_(blocks), jop_(jop) {
  const int norb = jop->nocc();
  Matrix mo1e = *jop->mo1e()->matrix();
  const btas::TensorView4<double> mo2e = btas::make_view(btas::CRange<4>(norb, norb, norb, norb), jop->mo2e()->storage());

  const int rasnorb = norb - (blocks->left_block()->norb() + blocks->right_block()->norb());
  const int lnorb = blocks->left_block()->norb();
  const int rnorb = blocks->right_block()->norb();

  { // make left jop to be used for left BlockOperators1
    const int leftsize = rasnorb+lnorb;
    auto left1e = make_shared<CSymMatrix>(mo1e.get_submatrix(0, 0, leftsize, leftsize));
    auto left2e = make_shared<Matrix>(leftsize*leftsize, leftsize*leftsize);
    auto low = {0, 0, 0, 0};
    auto up = {leftsize, leftsize, leftsize, leftsize};
    const btas::TensorView4<double> mo2eslice = btas::make_view(mo2e.range().slice(low, up), mo2e.storage());
    copy(mo2eslice.begin(), mo2eslice.end(), left2e->begin());
    auto leftjop = make_shared<DimerJop>(rasnorb, lnorb, left1e, left2e);
    left_ops_ = make_shared<BlockOperators1>(blocks->left_block(), leftjop);
  }

  { // make double block one to be used for intra BlockOperators1
    const int blocksize = lnorb + rnorb;
    auto block1e = make_shared<CSymMatrix>(mo1e.get_submatrix(rasnorb, rasnorb, blocksize, blocksize));
    auto block2e = make_shared<Matrix>(blocksize*blocksize, blocksize*blocksize);
    auto low = {rasnorb, rasnorb, rasnorb, rasnorb};
    auto up = {norb, norb, norb, norb};
    const btas::TensorView4<double> mo2eslice = btas::make_view(mo2e.range().slice(low, up), mo2e.storage());
    copy(mo2eslice.begin(), mo2eslice.end(), block2e->begin());
    auto intrajop = make_shared<DimerJop>(lnorb, rnorb, block1e, block2e);
    intra_ops_ = make_shared<BlockOperators1>(blocks->right_block(), intrajop);
  }

  { // make right jop to be used for right BlockOperators1
    const int rightsize = rasnorb + rnorb;
    auto right1efull = mo1e.get_submatrix(0, 0, rightsize, rightsize); // keeps the rasnorb x rasnorb terms in there
    right1efull->copy_block(0, rasnorb, rasnorb, rnorb, mo1e.get_submatrix(0, rasnorb+lnorb, rasnorb, rnorb));
    right1efull->copy_block(rasnorb, 0, rnorb, rasnorb, mo1e.get_submatrix(rasnorb+lnorb, 0, rnorb, rasnorb));
    right1efull->copy_block(rasnorb, rasnorb, rnorb, rnorb, mo1e.get_submatrix(rasnorb+lnorb, rasnorb+lnorb, rnorb, rnorb));
    auto right1e = make_shared<CSymMatrix>(right1efull);

    auto right2e = make_shared<Matrix>(rightsize*rightsize, rightsize*rightsize);
    btas::TensorView4<double> right2eView = btas::make_rwview(btas::CRange<4>(rightsize, rightsize, rightsize, rightsize), right2e->storage());
    for (int i = 0; i < rightsize; ++i) {
      const int ii = i < rasnorb ? i : i + lnorb;
      for (int j = 0; j < rightsize; ++j) {
        const int jj = j < rasnorb ? j : j + lnorb;
        for (int k = 0; k < rightsize; ++k) {
          const int kk = k < rasnorb ? k : k + lnorb;
          for (int l = 0; l < rightsize; ++l) {
            const int ll = l < rasnorb ? l : l + lnorb;
            right2eView(l, k, j, i) = mo2e(ll, kk, jj, ii);
          }
        }
      }
    }
    auto rightjop = make_shared<DimerJop>(rasnorb, rnorb, right1e, right2e);
    right_ops_ = make_shared<BlockOperators1>(blocks->right_block(), rightjop);
  }
}

// TODO: use sparsity
shared_ptr<BlockSparseMatrix> BlockOperators2::gamma_a(const BlockKey bk, int i) const {
  const int lnorb = blocks_->left_block()->norb();
  const BlockInfo binfo = blocks_->blockinfo(bk);
  const BlockKey target_bk(bk.nelea+1, bk.neleb);
  assert(blocks_->contains(target_bk));
  const BlockInfo tinfo = blocks_->blockinfo(target_bk);

  auto out = make_shared<Matrix>(tinfo.nstates, binfo.nstates, true);

  const vector<DMRG::BlockPair>& target_pairs = blocks_->blockpairs(target_bk);
  if (i < lnorb) { // i in L block
    for (auto& spair : blocks_->blockpairs(bk)) {
      BlockInfo rblock = spair.right;
      BlockInfo lsource = spair.left;

      auto iter = find_if(target_pairs.begin(), target_pairs.end(), [&rblock, &lsource] (DMRG::BlockPair bp)
        { return make_pair(BlockKey(lsource.nelea+1, lsource.neleb), rblock.key())==make_pair(bp.left.key(), bp.right.key()); } );
      if (iter != target_pairs.end()) {
        Matrix Runit(rblock.nstates, rblock.nstates); Runit.unit();
        Matrix Lgamma = *left_ops_->gamma_a_as_matrix(lsource.key(), i);

        out->add_block(1.0, iter->offset, spair.offset, iter->nstates(), spair.nstates(), kronecker_product(false, Runit, false, Lgamma));
      }
    }
  }
  else {
    for (auto& spair : blocks_->blockpairs(bk)) {
      BlockInfo lblock = spair.left;
      BlockInfo rsource = spair.right;

      const int phase = 1 - (((lblock.nelea + lblock.neleb)%2) << 1);

      auto iter = find_if(target_pairs.begin(), target_pairs.end(), [&lblock, &rsource] (DMRG::BlockPair bp)
        { return make_pair(lblock.key(), BlockKey(rsource.nelea+1, rsource.neleb))==make_pair(bp.left.key(), bp.right.key()); } );
      if (iter != target_pairs.end()) {
        Matrix Lunit(lblock.nstates, lblock.nstates); Lunit.unit();
        Matrix Rgamma = *right_ops_->gamma_a_as_matrix(rsource.key(), i - lnorb);

        out->add_block(phase, iter->offset, spair.offset, iter->nstates(), spair.nstates(), kronecker_product(false, Rgamma, false, Lunit));
      }
    }
  }

  return make_shared<BlockSparseMatrix>(out);
}

// TODO: use sparsity
shared_ptr<BlockSparseMatrix> BlockOperators2::gamma_b(const BlockKey bk, int i) const {
  const int lnorb = blocks_->left_block()->norb();
  const BlockInfo binfo = blocks_->blockinfo(bk);
  const BlockKey target_bk(bk.nelea, bk.neleb+1);
  assert(blocks_->contains(target_bk));
  const BlockInfo tinfo = blocks_->blockinfo(target_bk);

  auto out = make_shared<Matrix>(tinfo.nstates, binfo.nstates, true);

  const vector<DMRG::BlockPair>& target_pairs = blocks_->blockpairs(target_bk);
  if (i < lnorb) { // i in L block
    for (auto& spair : blocks_->blockpairs(bk)) {
      BlockInfo rblock = spair.right;
      BlockInfo lsource = spair.left;

      auto iter = find_if(target_pairs.begin(), target_pairs.end(), [&rblock, &lsource] (DMRG::BlockPair bp)
        { return make_pair(BlockKey(lsource.nelea, lsource.neleb+1), rblock.key())==make_pair(bp.left.key(), bp.right.key()); } );
      if (iter != target_pairs.end()) {
        Matrix Runit(rblock.nstates, rblock.nstates); Runit.unit();
        Matrix Lgamma = *left_ops_->gamma_b_as_matrix(lsource.key(), i);

        out->add_block(1.0, iter->offset, spair.offset, iter->nstates(), spair.nstates(), kronecker_product(false, Runit, false, Lgamma));
      }
    }
  }
  else {
    for (auto& spair : blocks_->blockpairs(bk)) {
      BlockInfo lblock = spair.left;
      BlockInfo rsource = spair.right;

      const int phase = 1 - (((lblock.nelea + lblock.neleb)%2) << 1);

      auto iter = find_if(target_pairs.begin(), target_pairs.end(), [&lblock, &rsource] (DMRG::BlockPair bp)
        { return make_pair(lblock.key(), BlockKey(rsource.nelea, rsource.neleb+1))==make_pair(bp.left.key(), bp.right.key()); } );
      if (iter != target_pairs.end()) {
        Matrix Lunit(lblock.nstates, lblock.nstates); Lunit.unit();
        Matrix Rgamma = *right_ops_->gamma_b_as_matrix(rsource.key(), i - lnorb);

        out->add_block(phase, iter->offset, spair.offset, iter->nstates(), spair.nstates(), kronecker_product(false, Rgamma, false, Lunit));
      }
    }
  }

  return make_shared<BlockSparseMatrix>(out);
}

double BlockOperators2::D_a(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const {
  const DMRG::BlockPair source_pair = *find_if(blocks_->blockpairs(bk).begin(), blocks_->blockpairs(bk).end(), [ketstate] (const DMRG::BlockPair bp)
    { return (ketstate >= bp.offset && ketstate < bp.offset+bp.nstates()); } );

  const BlockKey target_key(bk.nelea+1, bk.neleb);
  const DMRG::BlockPair target_pair = *find_if(blocks_->blockpairs(target_key).begin(), blocks_->blockpairs(target_key).end(), [brastate] (const DMRG::BlockPair bp)
    { return (brastate >= bp.offset && brastate < bp.offset+bp.nstates()); } );

  const int left_ketstate = (ketstate - source_pair.offset) % source_pair.left.nstates;
  const int right_ketstate = (ketstate - source_pair.offset) / source_pair.left.nstates;

  const int left_brastate = (brastate - target_pair.offset) % target_pair.left.nstates;
  const int right_brastate = (brastate - target_pair.offset) / target_pair.left.nstates;

  double out = 0.0;

  if (source_pair.left.key() == target_pair.left.key()) {
    if (left_ketstate==left_brastate)
      out += right_ops_->D_a(source_pair.right.key(), right_brastate, right_ketstate, i, j, k);
  }
  else if (source_pair.right.key() == target_pair.right.key()) {
    if (right_ketstate==right_brastate)
      out += left_ops_->D_a(source_pair.left.key(), left_brastate, left_ketstate, i, j, k);
  }
  return out;
}

double BlockOperators2::D_b(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const {
  const DMRG::BlockPair source_pair = *find_if(blocks_->blockpairs(bk).begin(), blocks_->blockpairs(bk).end(), [ketstate] (const DMRG::BlockPair bp)
    { return (ketstate >= bp.offset && ketstate < bp.offset+bp.nstates()); } );

  const BlockKey target_key(bk.nelea, bk.neleb+1);
  const DMRG::BlockPair target_pair = *find_if(blocks_->blockpairs(target_key).begin(), blocks_->blockpairs(target_key).end(), [brastate] (const DMRG::BlockPair bp)
    { return (brastate >= bp.offset && brastate < bp.offset+bp.nstates()); } );

  const int left_ketstate = (ketstate - source_pair.offset) % source_pair.left.nstates;
  const int right_ketstate = (ketstate - source_pair.offset) / source_pair.left.nstates;

  const int left_brastate = (brastate - target_pair.offset) % target_pair.left.nstates;
  const int right_brastate = (brastate - target_pair.offset) / target_pair.left.nstates;

  double out = 0.0;

  if (source_pair.left.key() == target_pair.left.key()) {
    if (left_ketstate==left_brastate)
      out += right_ops_->D_b(source_pair.right.key(), right_brastate, right_ketstate, i, j, k);
  }
  else if (source_pair.right.key() == target_pair.right.key()) {
    if (right_ketstate==right_brastate)
      out += left_ops_->D_b(source_pair.left.key(), left_brastate, left_ketstate, i, j, k);
  }
  return out;
}

// all of the offdiagonal blocks get an extra factor of two so I can call H = 0.5 * (H + H^t) at the end
shared_ptr<Matrix> BlockOperators2::ham(const BlockKey bk) const {
  const vector<DMRG::BlockPair>& block_pairs = blocks_->blockpairs(bk);

  BlockInfo binfo = blocks_->blockinfo(bk);

  auto out = make_shared<Matrix>(binfo.nstates, binfo.nstates);

  const int lnorb = blocks_->left_block()->norb();

  for (auto& source_pair : block_pairs) {
    { // diag parts
      // H_A (x) I_B + I_A (x) H_B
      shared_ptr<const Matrix> leftham = left_ops_->ham(source_pair.left);
      Matrix leftunit = *leftham->clone(); leftunit.unit();
      shared_ptr<const Matrix> rightham = right_ops_->ham(source_pair.right);
      Matrix rightunit = *rightham->clone(); rightunit.unit();

      Matrix diag = kronecker_product(false, rightunit, false, *leftham) + kronecker_product(false, *rightham, false, leftunit);

      // Q_aa (x) A^t   A
      const MatView Q_aa_view = intra_ops_->Q_aa_as_matview(source_pair.right.key());
      const MatView gamma_aa_view(btas::make_view(btas::CRange<2>(source_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb),
                               blocks_->left_block()->coupling({GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({source_pair.left.key(),source_pair.left.key()}).data->storage()), true);

      const MatView Q_bb_view = intra_ops_->Q_bb_as_matview(source_pair.right.key());
      const MatView gamma_bb_view(btas::make_view(btas::CRange<2>(source_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb),
                               blocks_->left_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({source_pair.left.key(),source_pair.left.key()}).data->storage()), true);

      Matrix ham_block = (gamma_aa_view ^ Q_aa_view) + (gamma_bb_view ^ Q_bb_view);

      sort_indices<0,2,1,3,1,1,1,1>(ham_block.data(), diag.data(), source_pair.left.nstates, source_pair.left.nstates, source_pair.right.nstates, source_pair.right.nstates);

      out->add_block(1.0, source_pair.offset, source_pair.offset, source_pair.nstates(), source_pair.nstates(), diag);
    }

    { // Q_ab (x) A^t   B
      const BlockKey right_target(source_pair.right.nelea-1, source_pair.right.neleb+1);
      const BlockKey left_target(source_pair.left.nelea+1, source_pair.left.neleb-1);

      auto iter = find_if(block_pairs.begin(), block_pairs.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );

      if (iter != block_pairs.end()) {
        DMRG::BlockPair target_pair = *iter;

        const MatView Q_ab_view = intra_ops_->Q_ab_as_matview(source_pair.right.key());
        const MatView gamma_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb),
                                 blocks_->left_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}).at({source_pair.left.key(),left_target}).data->storage()), true);

        // TODO maybe it would make more sense to reorder the block_ops in the first place?
        Matrix Qab(Q_ab_view);
        sort_indices<0,2,1,0,1,1,1>(Q_ab_view.data(), Qab.data(), Qab.ndim(), lnorb, lnorb);

        Matrix ham_block = gamma_view ^ Qab;

        Matrix tmp(target_pair.left.nstates*target_pair.right.nstates, source_pair.left.nstates*source_pair.right.nstates);
        sort_indices<1,2,0,3,0,1,1,1>(ham_block.data(), tmp.data(), source_pair.left.nstates, target_pair.left.nstates, target_pair.right.nstates, source_pair.right.nstates);

        out->add_block(2.0, target_pair.offset, source_pair.offset, tmp.ndim(), tmp.mdim(), tmp);
      }
    }

    { // P_aa (x) A^t A^t
      const BlockKey right_target(source_pair.right.nelea-2, source_pair.right.neleb);
      const BlockKey left_target(source_pair.left.nelea+2, source_pair.left.neleb);

      auto iter = find_if(block_pairs.begin(), block_pairs.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );

      if (iter != block_pairs.end()) {
        DMRG::BlockPair target_pair = *iter;

        const MatView P_aa_view = intra_ops_->P_aa_as_matview(source_pair.right.key());
        const MatView gamma_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb),
                                 blocks_->left_block()->coupling({GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}).at({source_pair.left.key(),left_target}).data->storage()), true);

        Matrix ham_block = gamma_view ^ P_aa_view;
        Matrix tmp(target_pair.nstates(), source_pair.nstates());

        sort_indices<1,2,0,3,0,1,-1,1>(ham_block.data(), tmp.data() , source_pair.left.nstates, target_pair.left.nstates, target_pair.right.nstates, source_pair.right.nstates);
        out->add_block(2.0, target_pair.offset, source_pair.offset, tmp.ndim(), tmp.mdim(), tmp);
      }
    }

    { // P_bb (x) B^t B^t
      const BlockKey right_target(source_pair.right.nelea, source_pair.right.neleb-2);
      const BlockKey left_target(source_pair.left.nelea, source_pair.left.neleb+2);

      auto iter = find_if(block_pairs.begin(), block_pairs.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );

      if (iter != block_pairs.end()) {
        DMRG::BlockPair target_pair = *iter;

        const MatView P_bb_view = intra_ops_->P_bb_as_matview(source_pair.right.key());
        const MatView gamma_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb),
                                 blocks_->left_block()->coupling({GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}).at({source_pair.left.key(),left_target}).data->storage()), true);

        Matrix ham_block = P_bb_view ^ gamma_view;
        Matrix tmp(target_pair.nstates(), source_pair.nstates());

        sort_indices<1,2,0,3,0,1,-1,1>(ham_block.data(), tmp.data() , source_pair.left.nstates, target_pair.left.nstates, target_pair.right.nstates, source_pair.right.nstates);
        out->add_block(2.0, target_pair.offset, source_pair.offset, tmp.ndim(), tmp.mdim(), tmp);
      }
    }

    { // P_ab (x) A^t B^t
      const BlockKey right_target(source_pair.right.nelea-1, source_pair.right.neleb-1);
      const BlockKey left_target(source_pair.left.nelea+1, source_pair.left.neleb+1);

      auto iter = find_if(block_pairs.begin(), block_pairs.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );

      if (iter != block_pairs.end()) {
        DMRG::BlockPair target_pair = *iter;

        const MatView P_ab_view = intra_ops_->P_ab_as_matview(source_pair.right.key());
        const MatView gamma_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb),
                                 blocks_->left_block()->coupling({GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}).at({source_pair.left.key(),left_target}).data->storage()), true);

        Matrix Pab(P_ab_view);
        sort_indices<0,2,1,0,1,1,1>(P_ab_view.data(), Pab.data(), Pab.ndim(), lnorb, lnorb);

        Matrix ham_block = gamma_view ^ Pab;

        Matrix tmp(target_pair.nstates(), source_pair.nstates());
        sort_indices<1,2,0,3,0,1,1,1>(ham_block.data(), tmp.data(), source_pair.left.nstates, target_pair.left.nstates, target_pair.right.nstates, source_pair.right.nstates);

        out->add_block(2.0, target_pair.offset, source_pair.offset, tmp.ndim(), tmp.mdim(), tmp);
      }
    }

    { // S_a (x)   A   and    D_a (x)   A^t A A + B^t B A
      const BlockKey left_target(source_pair.left.nelea-1, source_pair.left.neleb);
      const BlockKey right_target(source_pair.right.nelea+1, source_pair.right.neleb);
      auto iter = find_if(block_pairs.begin(), block_pairs.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );

      if (iter != block_pairs.end()) {
        DMRG::BlockPair target_pair = *iter;

        const MatView S_a_view = intra_ops_->S_a_as_matview(source_pair.right.key());
        const MatView gamma_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb),
                                 blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({source_pair.left.key(),left_target}).data->storage()), true);


        const MatView D_a_view = intra_ops_->D_a_as_matview(source_pair.right.key());
        shared_ptr<const btas::Tensor3<double>> gamma_aaa = blocks_->left_block()->coupling({GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({source_pair.left.key(),left_target}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_abb = blocks_->left_block()->coupling({GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({source_pair.left.key(),left_target}).data;
        const MatView gamma_aaa_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb*lnorb), gamma_aaa->storage()), true);
        const MatView gamma_abb_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb*lnorb), gamma_abb->storage()), true);

        Matrix gamma_akk(gamma_aaa_view + gamma_abb_view);

        Matrix Da(D_a_view);
        sort_indices<0,3,2,1,0,1,1,1>(D_a_view.data(), Da.data(), Da.ndim(), lnorb, lnorb, lnorb);

        Matrix ham_block = (gamma_view ^ S_a_view) + (gamma_akk ^ Da);

        Matrix tmp(target_pair.nstates(), source_pair.nstates());
        sort_indices<1,2,0,3,0,1,1,1>(ham_block.data(), tmp.data(), source_pair.left.nstates, target_pair.left.nstates, target_pair.right.nstates, source_pair.right.nstates);

        const double fac = -2.0 * static_cast<int>(1 - (((source_pair.left.nelea+source_pair.left.neleb)%2) << 1));
        out->add_block(fac, target_pair.offset, source_pair.offset, tmp.ndim(), tmp.mdim(), tmp);
      }
    }

    { // S_b (x)   B   and    D_b (x)   B^t B B + A^t A B
      const BlockKey left_target(source_pair.left.nelea, source_pair.left.neleb-1);
      const BlockKey right_target(source_pair.right.nelea, source_pair.right.neleb+1);
      auto iter = find_if(block_pairs.begin(), block_pairs.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );

      if (iter != block_pairs.end()) {
        DMRG::BlockPair target_pair = *iter;

        const MatView S_b_view = intra_ops_->S_b_as_matview(source_pair.right.key());
        const MatView gamma_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb),
                                 blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({source_pair.left.key(),left_target}).data->storage()), true);


        const MatView D_b_view = intra_ops_->D_b_as_matview(source_pair.right.key());
        shared_ptr<const btas::Tensor3<double>> gamma_bbb = blocks_->left_block()->coupling({GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({source_pair.left.key(),left_target}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_baa = blocks_->left_block()->coupling({GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({source_pair.left.key(),left_target}).data;
        const MatView gamma_bbb_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb*lnorb), gamma_bbb->storage()), true);
        const MatView gamma_baa_view(btas::make_view(btas::CRange<2>(target_pair.left.nstates*source_pair.left.nstates, lnorb*lnorb*lnorb), gamma_baa->storage()), true);

        Matrix gamma_bkk(gamma_bbb_view + gamma_baa_view);

        Matrix Db(D_b_view);
        sort_indices<0,3,2,1,0,1,1,1>(D_b_view.data(), Db.data(), Db.ndim(), lnorb, lnorb, lnorb);

        Matrix ham_block = (gamma_view ^ S_b_view) + (gamma_bkk ^ Db);

        Matrix tmp(target_pair.nstates(), source_pair.nstates());
        sort_indices<1,2,0,3,0,1,1,1>(ham_block.data(), tmp.data(), source_pair.left.nstates, target_pair.left.nstates, target_pair.right.nstates, source_pair.right.nstates);

        const double fac = -2.0 * static_cast<int>(1 - (((source_pair.left.nelea+source_pair.left.neleb)%2) << 1));
        out->add_block(fac, target_pair.offset, source_pair.offset, tmp.ndim(), tmp.mdim(), tmp);
      }
    }
  }

  // because symmetrize is actually very slow right now
  const size_t dim = out->ndim();
  double* data = out->data();
  for (size_t i = 0; i < dim; ++i)
    for (size_t j = 0; j < i; ++j)
      data[j + i*dim] = data[i + j*dim] = 0.5 * (data[j+i*dim] + data[i + j*dim]);

  return out;
}

/*-----------------------------------------------------------------------------------------------------
  All code below this point was generated by gen/blockops.py
-------------------------------------------------------------------------------------------------------*/

shared_ptr<BlockSparseMatrix> BlockOperators2::S_a(BlockKey bk, const int i) const {
  BlockKey target_bk(bk.nelea+1,bk.neleb);
  assert(blocks_->contains(target_bk));

  const vector<DMRG::BlockPair>& source_pvec = blocks_->blockpairs(bk);
  const vector<DMRG::BlockPair>& target_pvec = blocks_->blockpairs(target_bk);

  const BlockInfo source_info = blocks_->blockinfo(bk);
  const BlockInfo target_info = blocks_->blockinfo(target_bk);

  const int norb = jop_->nocc();
  const int lnorb = blocks_->left_block()->norb();
  const int rnorb = blocks_->right_block()->norb();
  const int loffset = norb - (lnorb + rnorb); // convenience variable for offset of left orbitals from zero
  const int roffset = norb - rnorb; // convenience variable for offset of right orbitals from zero

  const double* mo2e = jop_->mo2e()->data();

  const int max_L_M = max_element(blocks_->left_block()->blocks().begin(), blocks_->left_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const int max_R_M = max_element(blocks_->right_block()->blocks().begin(), blocks_->right_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const size_t max_L_intermediate = max_L_M * max_L_M * rnorb;
  const size_t max_R_intermediate = max_R_M * max_R_M * lnorb;
  const size_t max_intermediate = max(max_L_intermediate, max_R_intermediate);
  unique_ptr<double[]> scratch(new double[max_intermediate]);

  const size_t max_coulomb_size = lnorb < rnorb ? rnorb*rnorb*lnorb : lnorb*lnorb*rnorb;
  unique_ptr<double[]> coulomb(new double[max_coulomb_size]);
  map<pair<size_t, size_t>, shared_ptr<Matrix>> out;

  for (auto& spair : source_pvec) {
    // phase accumulated by moving an operator past the whole left ket block
    const int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%2) << 1);

    { // - 1.0 <L'| A^t B |L> (x) <R'| B^t |R>
      const BlockKey left_target(spair.left.nelea+1, spair.left.neleb-1);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb+1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({tpair.right.key(),spair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            for (int b = 0; b < lnorb; ++b) {
              coulomb[a + b*lnorb + p*lnorb*lnorb] = mo2e[(i)+(p+roffset)*norb+(b+loffset)*norb*norb+(a+loffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| B |L> (x) <R'| A^t B^t |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb-1);
      const BlockKey right_target(spair.right.nelea+1, spair.right.neleb+1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= lnorb*rnorb*rnorb);
        fill_n(scratch.get(), Rsize * lnorb, 0.0);
        for (int a = 0; a < lnorb; ++a) {
          for (int p = 0; p < rnorb; ++p) {
            for (int q = 0; q < rnorb; ++q) {
              coulomb[p + q*rnorb + a*rnorb*rnorb] = mo2e[(i)+(p+roffset)*norb+(q+roffset)*norb*norb+(a+loffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, 1.0, Rgamma->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
        for (int a = 0; a < lnorb; ++a) {
          const double* ldata = Lgamma->data() + Lsize * a;
          const double* rdata = scratch.get() + Rsize * a;
          kronecker_product(1.0, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| A |L> (x) <R'| A^t A^t |R>
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea+2, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= lnorb*rnorb*rnorb);
        fill_n(scratch.get(), Rsize * lnorb, 0.0);
        for (int a = 0; a < lnorb; ++a) {
          for (int p = 0; p < rnorb; ++p) {
            for (int q = 0; q < rnorb; ++q) {
              coulomb[p + q*rnorb + a*rnorb*rnorb] = mo2e[(i)+(p+roffset)*norb+(q+roffset)*norb*norb+(a+loffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, 1.0, Rgamma->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
        for (int a = 0; a < lnorb; ++a) {
          const double* ldata = Lgamma->data() + Lsize * a;
          const double* rdata = scratch.get() + Rsize * a;
          kronecker_product(1.0, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // S_a (x) I  + 1.0 <L'| A^t |L> (x) <R'| A^t A |R> - 1.0 <L'| A^t |L> (x) <R'| A^t A |R> + 1.0 <L'| A^t |L> (x) <R'| B^t B |R>
      const BlockKey left_target(spair.left.nelea+1, spair.left.neleb);
      const BlockKey right_target = spair.right.key();

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // S_a (x) I
        left_ops_->S_a_copy_to(scratch.get(), spair.left.key(), i);
        kronecker_product_I_B(1.0, spair.right.nstates, false, tpair.left.nstates, spair.left.nstates, scratch.get(), tpair.left.nstates, out_block->data(), out_block->ndim());

        { // ["1.0 <L'| A^t |L> (x) <R'| A^t A |R>", "1.0 <L'| A^t |L> (x) <R'| B^t B |R>"]
          shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({tpair.left.key(),spair.left.key()}).data;
          shared_ptr<const btas::Tensor3<double>> Rgamma1 = blocks_->right_block()->coupling({GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({tpair.right.key(),spair.right.key()}).data;
          shared_ptr<const btas::Tensor3<double>> Rgamma2 = blocks_->right_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({tpair.right.key(),spair.right.key()}).data;

          const int Lndim = Lgamma->extent(0);
          const int Lmdim = Lgamma->extent(1);
          const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

          const int Rndim = Rgamma1->extent(0);
          const int Rmdim = Rgamma1->extent(1);
          const int Rsize = Rgamma1->extent(0) * Rgamma1->extent(1);

          // fill in coulomb part
          assert(max_coulomb_size >= lnorb*rnorb*rnorb);
          fill_n(scratch.get(), Rsize * lnorb, 0.0);
          for (int a = 0; a < lnorb; ++a) {
            for (int p = 0; p < rnorb; ++p) {
              for (int q = 0; q < rnorb; ++q) {
                coulomb[p + q*rnorb + a*rnorb*rnorb] = (mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(q+roffset)*norb*norb*norb] - mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(q+roffset)*norb*norb*norb]);
              }
            }
          }

          dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, 1.0, Rgamma1->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
          for (int a = 0; a < lnorb; ++a) {
            for (int p = 0; p < rnorb; ++p) {
              for (int q = 0; q < rnorb; ++q) {
                coulomb[p + q*rnorb + a*rnorb*rnorb] = mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(q+roffset)*norb*norb*norb];
              }
            }
          }

          dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, 1.0, Rgamma2->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
          for (int a = 0; a < lnorb; ++a) {
            const double* ldata = Lgamma->data() + Lsize * a;
            const double* rdata = scratch.get() + Rsize * a;
            kronecker_product(1.0, false, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
          }
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| A^t A^t |L> (x) <R'| A |R>
      const BlockKey left_target(spair.left.nelea+2, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            for (int b = 0; b < lnorb; ++b) {
              coulomb[a + b*lnorb + p*lnorb*lnorb] = mo2e[(i)+(a+loffset)*norb+(b+loffset)*norb*norb+(p+roffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // I (x) S_a - 1.0 <L'| A^t A |L> (x) <R'| A^t |R> + 1.0 <L'| A^t A |L> (x) <R'| A^t |R> + 1.0 <L'| B^t B |L> (x) <R'| A^t |R>
      const BlockKey left_target = spair.left.key();
      const BlockKey right_target(spair.right.nelea+1, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // I (x) S_a
        right_ops_->S_a_copy_to(scratch.get(), spair.right.key(), i);
        kronecker_product_A_I(left_phase, false, tpair.right.nstates, spair.right.nstates, scratch.get(), tpair.right.nstates, spair.left.nstates, out_block->data(), out_block->ndim());

        { // ["-1.0 <L'| A^t A |L> (x) <R'| A^t |R>", "1.0 <L'| B^t B |L> (x) <R'| A^t |R>"]
          shared_ptr<const btas::Tensor3<double>> Lgamma1 = blocks_->left_block()->coupling({GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({tpair.left.key(),spair.left.key()}).data;
          shared_ptr<const btas::Tensor3<double>> Lgamma2 = blocks_->left_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({tpair.left.key(),spair.left.key()}).data;
          shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({tpair.right.key(),spair.right.key()}).data;

          const int Lndim = Lgamma1->extent(0);
          const int Lmdim = Lgamma1->extent(1);
          const int Lsize = Lgamma1->extent(0) * Lgamma1->extent(1);

          const int Rndim = Rgamma->extent(0);
          const int Rmdim = Rgamma->extent(1);
          const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

          // fill in coulomb part
          assert(max_coulomb_size >= rnorb*lnorb*lnorb);
          fill_n(scratch.get(), Lsize * rnorb, 0.0);
          for (int p = 0; p < rnorb; ++p) {
            for (int a = 0; a < lnorb; ++a) {
              for (int b = 0; b < lnorb; ++b) {
                coulomb[a + b*lnorb + p*lnorb*lnorb] = (mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(b+loffset)*norb*norb*norb] - mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(b+loffset)*norb*norb*norb]);
              }
            }
          }

          dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, -1.0, Lgamma1->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
          for (int p = 0; p < rnorb; ++p) {
            for (int a = 0; a < lnorb; ++a) {
              for (int b = 0; b < lnorb; ++b) {
                coulomb[a + b*lnorb + p*lnorb*lnorb] = mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(b+loffset)*norb*norb*norb];
              }
            }
          }

          dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, 1.0, Lgamma2->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
          for (int p = 0; p < rnorb; ++p) {
            const double* ldata = scratch.get() + Lsize * p;
            const double* rdata = Rgamma->data() + Rsize * p;
            kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
          }
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // - 1.0 <L'| B^t |L> (x) <R'| A^t B |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb+1);
      const BlockKey right_target(spair.right.nelea+1, spair.right.neleb-1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({tpair.left.key(),spair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= lnorb*rnorb*rnorb);
        fill_n(scratch.get(), Rsize * lnorb, 0.0);
        for (int a = 0; a < lnorb; ++a) {
          for (int p = 0; p < rnorb; ++p) {
            for (int q = 0; q < rnorb; ++q) {
              coulomb[p + q*rnorb + a*rnorb*rnorb] = mo2e[(i)+(a+loffset)*norb+(q+roffset)*norb*norb+(p+roffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, -1.0, Rgamma->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
        for (int a = 0; a < lnorb; ++a) {
          const double* ldata = Lgamma->data() + Lsize * a;
          const double* rdata = scratch.get() + Rsize * a;
          kronecker_product(1.0, true, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| A^t B^t |L> (x) <R'| B |R>
      const BlockKey left_target(spair.left.nelea+1, spair.left.neleb+1);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb-1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            for (int b = 0; b < lnorb; ++b) {
              coulomb[a + b*lnorb + p*lnorb*lnorb] = mo2e[(i)+(a+loffset)*norb+(b+loffset)*norb*norb+(p+roffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }
  }

  return make_shared<BlockSparseMatrix>(target_info.nstates, source_info.nstates, out);
}


shared_ptr<BlockSparseMatrix> BlockOperators2::S_b(BlockKey bk, const int i) const {
  BlockKey target_bk(bk.nelea,bk.neleb+1);
  assert(blocks_->contains(target_bk));

  const vector<DMRG::BlockPair>& source_pvec = blocks_->blockpairs(bk);
  const vector<DMRG::BlockPair>& target_pvec = blocks_->blockpairs(target_bk);

  const BlockInfo source_info = blocks_->blockinfo(bk);
  const BlockInfo target_info = blocks_->blockinfo(target_bk);

  const int norb = jop_->nocc();
  const int lnorb = blocks_->left_block()->norb();
  const int rnorb = blocks_->right_block()->norb();
  const int loffset = norb - (lnorb + rnorb); // convenience variable for offset of left orbitals from zero
  const int roffset = norb - rnorb; // convenience variable for offset of right orbitals from zero

  const double* mo2e = jop_->mo2e()->data();

  const int max_L_M = max_element(blocks_->left_block()->blocks().begin(), blocks_->left_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const int max_R_M = max_element(blocks_->right_block()->blocks().begin(), blocks_->right_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const size_t max_L_intermediate = max_L_M * max_L_M * rnorb;
  const size_t max_R_intermediate = max_R_M * max_R_M * lnorb;
  const size_t max_intermediate = max(max_L_intermediate, max_R_intermediate);
  unique_ptr<double[]> scratch(new double[max_intermediate]);

  const size_t max_coulomb_size = lnorb < rnorb ? rnorb*rnorb*lnorb : lnorb*lnorb*rnorb;
  unique_ptr<double[]> coulomb(new double[max_coulomb_size]);
  map<pair<size_t, size_t>, shared_ptr<Matrix>> out;

  for (auto& spair : source_pvec) {
    // phase accumulated by moving an operator past the whole left ket block
    const int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%2) << 1);

    { // - 1.0 <L'| B^t A |L> (x) <R'| A^t |R>
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb+1);
      const BlockKey right_target(spair.right.nelea+1, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}).at({tpair.left.key(),spair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({tpair.right.key(),spair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            for (int b = 0; b < lnorb; ++b) {
              coulomb[a + b*lnorb + p*lnorb*lnorb] = mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(b+loffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| B |L> (x) <R'| B^t B^t |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb-1);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb+2);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= lnorb*rnorb*rnorb);
        fill_n(scratch.get(), Rsize * lnorb, 0.0);
        for (int a = 0; a < lnorb; ++a) {
          for (int p = 0; p < rnorb; ++p) {
            for (int q = 0; q < rnorb; ++q) {
              coulomb[p + q*rnorb + a*rnorb*rnorb] = mo2e[(i)+(p+roffset)*norb+(q+roffset)*norb*norb+(a+loffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, 1.0, Rgamma->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
        for (int a = 0; a < lnorb; ++a) {
          const double* ldata = Lgamma->data() + Lsize * a;
          const double* rdata = scratch.get() + Rsize * a;
          kronecker_product(1.0, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // - 1.0 <L'| B^t A^t |L> (x) <R'| A |R>
      const BlockKey left_target(spair.left.nelea+1, spair.left.neleb+1);
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            for (int b = 0; b < lnorb; ++b) {
              coulomb[a + b*lnorb + p*lnorb*lnorb] = mo2e[(i)+(b+loffset)*norb+(a+loffset)*norb*norb+(p+roffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // S_b (x) I  + 1.0 <L'| B^t |L> (x) <R'| B^t B |R> - 1.0 <L'| B^t |L> (x) <R'| B^t B |R> + 1.0 <L'| B^t |L> (x) <R'| A^t A |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb+1);
      const BlockKey right_target = spair.right.key();

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // S_b (x) I
        left_ops_->S_b_copy_to(scratch.get(), spair.left.key(), i);
        kronecker_product_I_B(1.0, spair.right.nstates, false, tpair.left.nstates, spair.left.nstates, scratch.get(), tpair.left.nstates, out_block->data(), out_block->ndim());

        { // ["1.0 <L'| B^t |L> (x) <R'| B^t B |R>", "1.0 <L'| B^t |L> (x) <R'| A^t A |R>"]
          shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({tpair.left.key(),spair.left.key()}).data;
          shared_ptr<const btas::Tensor3<double>> Rgamma1 = blocks_->right_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({tpair.right.key(),spair.right.key()}).data;
          shared_ptr<const btas::Tensor3<double>> Rgamma2 = blocks_->right_block()->coupling({GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({tpair.right.key(),spair.right.key()}).data;

          const int Lndim = Lgamma->extent(0);
          const int Lmdim = Lgamma->extent(1);
          const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

          const int Rndim = Rgamma1->extent(0);
          const int Rmdim = Rgamma1->extent(1);
          const int Rsize = Rgamma1->extent(0) * Rgamma1->extent(1);

          // fill in coulomb part
          assert(max_coulomb_size >= lnorb*rnorb*rnorb);
          fill_n(scratch.get(), Rsize * lnorb, 0.0);
          for (int a = 0; a < lnorb; ++a) {
            for (int p = 0; p < rnorb; ++p) {
              for (int q = 0; q < rnorb; ++q) {
                coulomb[p + q*rnorb + a*rnorb*rnorb] = (mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(q+roffset)*norb*norb*norb] - mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(q+roffset)*norb*norb*norb]);
              }
            }
          }

          dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, 1.0, Rgamma1->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
          for (int a = 0; a < lnorb; ++a) {
            for (int p = 0; p < rnorb; ++p) {
              for (int q = 0; q < rnorb; ++q) {
                coulomb[p + q*rnorb + a*rnorb*rnorb] = mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(q+roffset)*norb*norb*norb];
              }
            }
          }

          dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, 1.0, Rgamma2->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
          for (int a = 0; a < lnorb; ++a) {
            const double* ldata = Lgamma->data() + Lsize * a;
            const double* rdata = scratch.get() + Rsize * a;
            kronecker_product(1.0, false, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
          }
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // - 1.0 <L'| A^t |L> (x) <R'| B^t A |R>
      const BlockKey left_target(spair.left.nelea+1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb+1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({tpair.left.key(),spair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}).at({tpair.right.key(),spair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= lnorb*rnorb*rnorb);
        fill_n(scratch.get(), Rsize * lnorb, 0.0);
        for (int a = 0; a < lnorb; ++a) {
          for (int p = 0; p < rnorb; ++p) {
            for (int q = 0; q < rnorb; ++q) {
              coulomb[p + q*rnorb + a*rnorb*rnorb] = mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(q+roffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, -1.0, Rgamma->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
        for (int a = 0; a < lnorb; ++a) {
          const double* ldata = Lgamma->data() + Lsize * a;
          const double* rdata = scratch.get() + Rsize * a;
          kronecker_product(1.0, false, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // - 1.0 <L'| A |L> (x) <R'| B^t A^t |R>
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea+1, spair.right.neleb+1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= lnorb*rnorb*rnorb);
        fill_n(scratch.get(), Rsize * lnorb, 0.0);
        for (int a = 0; a < lnorb; ++a) {
          for (int p = 0; p < rnorb; ++p) {
            for (int q = 0; q < rnorb; ++q) {
              coulomb[p + q*rnorb + a*rnorb*rnorb] = mo2e[(i)+(q+roffset)*norb+(p+roffset)*norb*norb+(a+loffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Rsize, lnorb, rnorb*rnorb, -1.0, Rgamma->data(), Rsize, coulomb.get(), rnorb*rnorb, 1.0, scratch.get(), Rsize);
        for (int a = 0; a < lnorb; ++a) {
          const double* ldata = Lgamma->data() + Lsize * a;
          const double* rdata = scratch.get() + Rsize * a;
          kronecker_product(1.0, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // I (x) S_b - 1.0 <L'| B^t B |L> (x) <R'| B^t |R> + 1.0 <L'| B^t B |L> (x) <R'| B^t |R> + 1.0 <L'| A^t A |L> (x) <R'| B^t |R>
      const BlockKey left_target = spair.left.key();
      const BlockKey right_target(spair.right.nelea, spair.right.neleb+1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // I (x) S_b
        right_ops_->S_b_copy_to(scratch.get(), spair.right.key(), i);
        kronecker_product_A_I(left_phase, false, tpair.right.nstates, spair.right.nstates, scratch.get(), tpair.right.nstates, spair.left.nstates, out_block->data(), out_block->ndim());

        { // ["-1.0 <L'| B^t B |L> (x) <R'| B^t |R>", "1.0 <L'| A^t A |L> (x) <R'| B^t |R>"]
          shared_ptr<const btas::Tensor3<double>> Lgamma1 = blocks_->left_block()->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({tpair.left.key(),spair.left.key()}).data;
          shared_ptr<const btas::Tensor3<double>> Lgamma2 = blocks_->left_block()->coupling({GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({tpair.left.key(),spair.left.key()}).data;
          shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({tpair.right.key(),spair.right.key()}).data;

          const int Lndim = Lgamma1->extent(0);
          const int Lmdim = Lgamma1->extent(1);
          const int Lsize = Lgamma1->extent(0) * Lgamma1->extent(1);

          const int Rndim = Rgamma->extent(0);
          const int Rmdim = Rgamma->extent(1);
          const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

          // fill in coulomb part
          assert(max_coulomb_size >= rnorb*lnorb*lnorb);
          fill_n(scratch.get(), Lsize * rnorb, 0.0);
          for (int p = 0; p < rnorb; ++p) {
            for (int a = 0; a < lnorb; ++a) {
              for (int b = 0; b < lnorb; ++b) {
                coulomb[a + b*lnorb + p*lnorb*lnorb] = (mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(b+loffset)*norb*norb*norb] - mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(b+loffset)*norb*norb*norb]);
              }
            }
          }

          dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, -1.0, Lgamma1->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
          for (int p = 0; p < rnorb; ++p) {
            for (int a = 0; a < lnorb; ++a) {
              for (int b = 0; b < lnorb; ++b) {
                coulomb[a + b*lnorb + p*lnorb*lnorb] = mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(b+loffset)*norb*norb*norb];
              }
            }
          }

          dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, 1.0, Lgamma2->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
          for (int p = 0; p < rnorb; ++p) {
            const double* ldata = scratch.get() + Lsize * p;
            const double* rdata = Rgamma->data() + Rsize * p;
            kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
          }
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| B^t B^t |L> (x) <R'| B |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb+2);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb-1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            for (int b = 0; b < lnorb; ++b) {
              coulomb[a + b*lnorb + p*lnorb*lnorb] = mo2e[(i)+(a+loffset)*norb+(b+loffset)*norb*norb+(p+roffset)*norb*norb*norb];
            }
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb*lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb*lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }
  }

  return make_shared<BlockSparseMatrix>(target_info.nstates, source_info.nstates, out);
}


shared_ptr<BlockSparseMatrix> BlockOperators2::Q_aa(BlockKey bk, const int i, const int j) const {
  const vector<DMRG::BlockPair>& pvec = blocks_->blockpairs(bk);

  const BlockInfo binfo = blocks_->blockinfo(bk);

  const int norb = jop_->nocc();
  const int lnorb = blocks_->left_block()->norb();
  const int rnorb = blocks_->right_block()->norb();
  const int loffset = norb - (lnorb + rnorb); // convenience variable for offset of left orbitals from zero
  const int roffset = norb - rnorb; // convenience variable for offset of right orbitals from zero

  const double* mo2e = jop_->mo2e()->data();

  const int max_L_M = max_element(blocks_->left_block()->blocks().begin(), blocks_->left_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const int max_R_M = max_element(blocks_->right_block()->blocks().begin(), blocks_->right_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const size_t max_L_intermediate = max_L_M * max_L_M * rnorb;
  const size_t max_R_intermediate = max_R_M * max_R_M * lnorb;
  const size_t max_intermediate = max(max_L_intermediate, max_R_intermediate);
  unique_ptr<double[]> scratch(new double[max_intermediate]);

  const size_t max_coulomb_size = rnorb*lnorb; // lnorb < rnorb ? rnorb*lnorb : lnorb*rnorb;
  unique_ptr<double[]> coulomb(new double[max_coulomb_size]);
  map<pair<size_t, size_t>, shared_ptr<Matrix>> out;

  for (auto& spair : pvec) {
    // phase accumulated by moving an operator past the whole left ket block
    const int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%2) << 1);

    { // I (x) Q_aa +  Q_aa (x) I
      auto out_block = make_shared<Matrix>(spair.nstates(), spair.nstates(), true);

      // I (x) Q_aa
      right_ops_->Q_aa_copy_to(scratch.get(), spair.right.key(), i, j);
      kronecker_product_A_I(1.0, false, spair.right.nstates, spair.right.nstates, scratch.get(), spair.right.nstates, spair.left.nstates, out_block->data(), out_block->ndim());

      // Q_aa (x) I
      left_ops_->Q_aa_copy_to(scratch.get(), spair.left.key(), i, j);
      kronecker_product_I_B(1.0, spair.right.nstates, false, spair.left.nstates, spair.left.nstates, scratch.get(), spair.left.nstates, out_block->data(), out_block->ndim());

      // add to map if large enough
      if (out_block->rms() > thresh_)
        out.emplace(make_pair<size_t, size_t>(spair.offset, spair.offset), out_block);
    }

    { // - 1.0 <L'| B |L> (x) <R'| B^t |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb-1);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb+1);

      auto iter = find_if(pvec.begin(), pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({tpair.right.key(),spair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = mo2e[(p+roffset)+(i)*norb+(a+loffset)*norb*norb+(j)*norb*norb*norb];
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // - 1.0 <L'| A |L> (x) <R'| A^t |R>
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea+1, spair.right.neleb);

      auto iter = find_if(pvec.begin(), pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({tpair.right.key(),spair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = (mo2e[(p+roffset)+(i)*norb+(a+loffset)*norb*norb+(j)*norb*norb*norb] - mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(j)*norb*norb*norb]);
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| B^t |L> (x) <R'| B |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb+1);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb-1);

      auto iter = find_if(pvec.begin(), pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({tpair.left.key(),spair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = mo2e[(a+loffset)+(i)*norb+(p+roffset)*norb*norb+(j)*norb*norb*norb];
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| A^t |L> (x) <R'| A |R>
      const BlockKey left_target(spair.left.nelea+1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb);

      auto iter = find_if(pvec.begin(), pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({tpair.left.key(),spair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = (mo2e[(a+loffset)+(i)*norb+(p+roffset)*norb*norb+(j)*norb*norb*norb] - mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(j)*norb*norb*norb]);
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }
  }

  return make_shared<BlockSparseMatrix>(binfo.nstates, binfo.nstates, out);
}


shared_ptr<BlockSparseMatrix> BlockOperators2::Q_bb(BlockKey bk, const int i, const int j) const {
  const vector<DMRG::BlockPair>& pvec = blocks_->blockpairs(bk);

  const BlockInfo binfo = blocks_->blockinfo(bk);

  const int norb = jop_->nocc();
  const int lnorb = blocks_->left_block()->norb();
  const int rnorb = blocks_->right_block()->norb();
  const int loffset = norb - (lnorb + rnorb); // convenience variable for offset of left orbitals from zero
  const int roffset = norb - rnorb; // convenience variable for offset of right orbitals from zero

  const double* mo2e = jop_->mo2e()->data();

  const int max_L_M = max_element(blocks_->left_block()->blocks().begin(), blocks_->left_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const int max_R_M = max_element(blocks_->right_block()->blocks().begin(), blocks_->right_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const size_t max_L_intermediate = max_L_M * max_L_M * rnorb;
  const size_t max_R_intermediate = max_R_M * max_R_M * lnorb;
  const size_t max_intermediate = max(max_L_intermediate, max_R_intermediate);
  unique_ptr<double[]> scratch(new double[max_intermediate]);

  const size_t max_coulomb_size = rnorb*lnorb; // lnorb < rnorb ? rnorb*lnorb : lnorb*rnorb;
  unique_ptr<double[]> coulomb(new double[max_coulomb_size]);
  map<pair<size_t, size_t>, shared_ptr<Matrix>> out;

  for (auto& spair : pvec) {
    // phase accumulated by moving an operator past the whole left ket block
    const int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%2) << 1);

    { // I (x) Q_bb +  Q_bb (x) I
      auto out_block = make_shared<Matrix>(spair.nstates(), spair.nstates(), true);

      // I (x) Q_bb
      right_ops_->Q_bb_copy_to(scratch.get(), spair.right.key(), i, j);
      kronecker_product_A_I(1.0, false, spair.right.nstates, spair.right.nstates, scratch.get(), spair.right.nstates, spair.left.nstates, out_block->data(), out_block->ndim());

      // Q_bb (x) I
      left_ops_->Q_bb_copy_to(scratch.get(), spair.left.key(), i, j);
      kronecker_product_I_B(1.0, spair.right.nstates, false, spair.left.nstates, spair.left.nstates, scratch.get(), spair.left.nstates, out_block->data(), out_block->ndim());

      // add to map if large enough
      if (out_block->rms() > thresh_)
        out.emplace(make_pair<size_t, size_t>(spair.offset, spair.offset), out_block);
    }

    { // - 1.0 <L'| B |L> (x) <R'| B^t |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb-1);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb+1);

      auto iter = find_if(pvec.begin(), pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({tpair.right.key(),spair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = (mo2e[(p+roffset)+(i)*norb+(a+loffset)*norb*norb+(j)*norb*norb*norb] - mo2e[(i)+(p+roffset)*norb+(a+loffset)*norb*norb+(j)*norb*norb*norb]);
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // - 1.0 <L'| A |L> (x) <R'| A^t |R>
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea+1, spair.right.neleb);

      auto iter = find_if(pvec.begin(), pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({tpair.right.key(),spair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = mo2e[(p+roffset)+(i)*norb+(a+loffset)*norb*norb+(j)*norb*norb*norb];
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| A^t |L> (x) <R'| A |R>
      const BlockKey left_target(spair.left.nelea+1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb);

      auto iter = find_if(pvec.begin(), pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({tpair.left.key(),spair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = mo2e[(a+loffset)+(i)*norb+(p+roffset)*norb*norb+(j)*norb*norb*norb];
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| B^t |L> (x) <R'| B |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb+1);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb-1);

      auto iter = find_if(pvec.begin(), pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({tpair.left.key(),spair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = (mo2e[(a+loffset)+(i)*norb+(p+roffset)*norb*norb+(j)*norb*norb*norb] - mo2e[(i)+(a+loffset)*norb+(p+roffset)*norb*norb+(j)*norb*norb*norb]);
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }
  }

  return make_shared<BlockSparseMatrix>(binfo.nstates, binfo.nstates, out);
}


shared_ptr<BlockSparseMatrix> BlockOperators2::Q_ab(BlockKey bk, const int i, const int j) const {
  BlockKey target_bk(bk.nelea-1,bk.neleb+1);
  assert(blocks_->contains(target_bk));

  const vector<DMRG::BlockPair>& source_pvec = blocks_->blockpairs(bk);
  const vector<DMRG::BlockPair>& target_pvec = blocks_->blockpairs(target_bk);

  const BlockInfo source_info = blocks_->blockinfo(bk);
  const BlockInfo target_info = blocks_->blockinfo(target_bk);

  const int norb = jop_->nocc();
  const int lnorb = blocks_->left_block()->norb();
  const int rnorb = blocks_->right_block()->norb();
  const int loffset = norb - (lnorb + rnorb); // convenience variable for offset of left orbitals from zero
  const int roffset = norb - rnorb; // convenience variable for offset of right orbitals from zero

  const double* mo2e = jop_->mo2e()->data();

  const int max_L_M = max_element(blocks_->left_block()->blocks().begin(), blocks_->left_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const int max_R_M = max_element(blocks_->right_block()->blocks().begin(), blocks_->right_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const size_t max_L_intermediate = max_L_M * max_L_M * rnorb;
  const size_t max_R_intermediate = max_R_M * max_R_M * lnorb;
  const size_t max_intermediate = max(max_L_intermediate, max_R_intermediate);
  unique_ptr<double[]> scratch(new double[max_intermediate]);

  const size_t max_coulomb_size = rnorb*lnorb; // lnorb < rnorb ? rnorb*lnorb : lnorb*rnorb;
  unique_ptr<double[]> coulomb(new double[max_coulomb_size]);
  map<pair<size_t, size_t>, shared_ptr<Matrix>> out;

  for (auto& spair : source_pvec) {
    // phase accumulated by moving an operator past the whole left ket block
    const int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%2) << 1);

    { // I (x) Q_ab
      const BlockKey left_target = spair.left.key();
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb+1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // I (x) Q_ab
        right_ops_->Q_ab_copy_to(scratch.get(), spair.right.key(), i, j);
        kronecker_product_A_I(1.0, false, tpair.right.nstates, spair.right.nstates, scratch.get(), tpair.right.nstates, spair.left.nstates, out_block->data(), out_block->ndim());

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| A |L> (x) <R'| B^t |R>
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb+1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({tpair.right.key(),spair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = mo2e[(p+roffset)+(a+loffset)*norb+(j)*norb*norb+(i)*norb*norb*norb];
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, false, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // Q_ab (x) I
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb+1);
      const BlockKey right_target = spair.right.key();

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // Q_ab (x) I
        left_ops_->Q_ab_copy_to(scratch.get(), spair.left.key(), i, j);
        kronecker_product_I_B(1.0, spair.right.nstates, false, tpair.left.nstates, spair.left.nstates, scratch.get(), tpair.left.nstates, out_block->data(), out_block->ndim());

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // - 1.0 <L'| B^t |L> (x) <R'| A |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb+1);
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({tpair.left.key(),spair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = mo2e[(a+loffset)+(p+roffset)*norb+(j)*norb*norb+(i)*norb*norb*norb];
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, false, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }
  }

  return make_shared<BlockSparseMatrix>(target_info.nstates, source_info.nstates, out);
}


shared_ptr<BlockSparseMatrix> BlockOperators2::P_aa(BlockKey bk, const int i, const int j) const {
  BlockKey target_bk(bk.nelea-2,bk.neleb);
  assert(blocks_->contains(target_bk));

  const vector<DMRG::BlockPair>& source_pvec = blocks_->blockpairs(bk);
  const vector<DMRG::BlockPair>& target_pvec = blocks_->blockpairs(target_bk);

  const BlockInfo source_info = blocks_->blockinfo(bk);
  const BlockInfo target_info = blocks_->blockinfo(target_bk);

  const int norb = jop_->nocc();
  const int lnorb = blocks_->left_block()->norb();
  const int rnorb = blocks_->right_block()->norb();
  const int loffset = norb - (lnorb + rnorb); // convenience variable for offset of left orbitals from zero
  const int roffset = norb - rnorb; // convenience variable for offset of right orbitals from zero

  const double* mo2e = jop_->mo2e()->data();

  const int max_L_M = max_element(blocks_->left_block()->blocks().begin(), blocks_->left_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const int max_R_M = max_element(blocks_->right_block()->blocks().begin(), blocks_->right_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const size_t max_L_intermediate = max_L_M * max_L_M * rnorb;
  const size_t max_R_intermediate = max_R_M * max_R_M * lnorb;
  const size_t max_intermediate = max(max_L_intermediate, max_R_intermediate);
  unique_ptr<double[]> scratch(new double[max_intermediate]);

  const size_t max_coulomb_size = rnorb*lnorb; // lnorb < rnorb ? rnorb*lnorb : lnorb*rnorb;
  unique_ptr<double[]> coulomb(new double[max_coulomb_size]);
  map<pair<size_t, size_t>, shared_ptr<Matrix>> out;

  for (auto& spair : source_pvec) {
    // phase accumulated by moving an operator past the whole left ket block
    const int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%2) << 1);

    { // + 0.5 <L'| A |L> (x) <R'| A |R> - 0.5 <L'| A |L> (x) <R'| A |R>
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = (mo2e[(a+loffset)+(p+roffset)*norb+(j)*norb*norb+(i)*norb*norb*norb] - mo2e[(a+loffset)+(p+roffset)*norb+(i)*norb*norb+(j)*norb*norb*norb]);
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, 0.5, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // I (x) P_aa
      const BlockKey left_target = spair.left.key();
      const BlockKey right_target(spair.right.nelea-2, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // I (x) P_aa
        right_ops_->P_aa_copy_to(scratch.get(), spair.right.key(), i, j);
        kronecker_product_A_I(1.0, false, tpair.right.nstates, spair.right.nstates, scratch.get(), tpair.right.nstates, spair.left.nstates, out_block->data(), out_block->ndim());

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // P_aa (x) I
      const BlockKey left_target(spair.left.nelea-2, spair.left.neleb);
      const BlockKey right_target = spair.right.key();

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // P_aa (x) I
        left_ops_->P_aa_copy_to(scratch.get(), spair.left.key(), i, j);
        kronecker_product_I_B(1.0, spair.right.nstates, false, tpair.left.nstates, spair.left.nstates, scratch.get(), tpair.left.nstates, out_block->data(), out_block->ndim());

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }
  }

  return make_shared<BlockSparseMatrix>(target_info.nstates, source_info.nstates, out);
}


shared_ptr<BlockSparseMatrix> BlockOperators2::P_bb(BlockKey bk, const int i, const int j) const {
  BlockKey target_bk(bk.nelea,bk.neleb-2);
  assert(blocks_->contains(target_bk));

  const vector<DMRG::BlockPair>& source_pvec = blocks_->blockpairs(bk);
  const vector<DMRG::BlockPair>& target_pvec = blocks_->blockpairs(target_bk);

  const BlockInfo source_info = blocks_->blockinfo(bk);
  const BlockInfo target_info = blocks_->blockinfo(target_bk);

  const int norb = jop_->nocc();
  const int lnorb = blocks_->left_block()->norb();
  const int rnorb = blocks_->right_block()->norb();
  const int loffset = norb - (lnorb + rnorb); // convenience variable for offset of left orbitals from zero
  const int roffset = norb - rnorb; // convenience variable for offset of right orbitals from zero

  const double* mo2e = jop_->mo2e()->data();

  const int max_L_M = max_element(blocks_->left_block()->blocks().begin(), blocks_->left_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const int max_R_M = max_element(blocks_->right_block()->blocks().begin(), blocks_->right_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const size_t max_L_intermediate = max_L_M * max_L_M * rnorb;
  const size_t max_R_intermediate = max_R_M * max_R_M * lnorb;
  const size_t max_intermediate = max(max_L_intermediate, max_R_intermediate);
  unique_ptr<double[]> scratch(new double[max_intermediate]);

  const size_t max_coulomb_size = rnorb*lnorb; // lnorb < rnorb ? rnorb*lnorb : lnorb*rnorb;
  unique_ptr<double[]> coulomb(new double[max_coulomb_size]);
  map<pair<size_t, size_t>, shared_ptr<Matrix>> out;

  for (auto& spair : source_pvec) {
    // phase accumulated by moving an operator past the whole left ket block
    const int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%2) << 1);

    { // + 0.5 <L'| B |L> (x) <R'| B |R> - 0.5 <L'| B |L> (x) <R'| B |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb-1);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb-1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = (mo2e[(a+loffset)+(p+roffset)*norb+(j)*norb*norb+(i)*norb*norb*norb] - mo2e[(a+loffset)+(p+roffset)*norb+(i)*norb*norb+(j)*norb*norb*norb]);
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, 0.5, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // I (x) P_bb
      const BlockKey left_target = spair.left.key();
      const BlockKey right_target(spair.right.nelea, spair.right.neleb-2);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // I (x) P_bb
        right_ops_->P_bb_copy_to(scratch.get(), spair.right.key(), i, j);
        kronecker_product_A_I(1.0, false, tpair.right.nstates, spair.right.nstates, scratch.get(), tpair.right.nstates, spair.left.nstates, out_block->data(), out_block->ndim());

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // P_bb (x) I
      const BlockKey left_target(spair.left.nelea, spair.left.neleb-2);
      const BlockKey right_target = spair.right.key();

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // P_bb (x) I
        left_ops_->P_bb_copy_to(scratch.get(), spair.left.key(), i, j);
        kronecker_product_I_B(1.0, spair.right.nstates, false, tpair.left.nstates, spair.left.nstates, scratch.get(), tpair.left.nstates, out_block->data(), out_block->ndim());

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }
  }

  return make_shared<BlockSparseMatrix>(target_info.nstates, source_info.nstates, out);
}


shared_ptr<BlockSparseMatrix> BlockOperators2::P_ab(BlockKey bk, const int i, const int j) const {
  BlockKey target_bk(bk.nelea-1,bk.neleb-1);
  assert(blocks_->contains(target_bk));

  const vector<DMRG::BlockPair>& source_pvec = blocks_->blockpairs(bk);
  const vector<DMRG::BlockPair>& target_pvec = blocks_->blockpairs(target_bk);

  const BlockInfo source_info = blocks_->blockinfo(bk);
  const BlockInfo target_info = blocks_->blockinfo(target_bk);

  const int norb = jop_->nocc();
  const int lnorb = blocks_->left_block()->norb();
  const int rnorb = blocks_->right_block()->norb();
  const int loffset = norb - (lnorb + rnorb); // convenience variable for offset of left orbitals from zero
  const int roffset = norb - rnorb; // convenience variable for offset of right orbitals from zero

  const double* mo2e = jop_->mo2e()->data();

  const int max_L_M = max_element(blocks_->left_block()->blocks().begin(), blocks_->left_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const int max_R_M = max_element(blocks_->right_block()->blocks().begin(), blocks_->right_block()->blocks().end(),
                                   [] (const BlockInfo& a, const BlockInfo& b)
                                     { return a.nstates < b.nstates; })->nstates;
  const size_t max_L_intermediate = max_L_M * max_L_M * rnorb;
  const size_t max_R_intermediate = max_R_M * max_R_M * lnorb;
  const size_t max_intermediate = max(max_L_intermediate, max_R_intermediate);
  unique_ptr<double[]> scratch(new double[max_intermediate]);

  const size_t max_coulomb_size = rnorb*lnorb; // lnorb < rnorb ? rnorb*lnorb : lnorb*rnorb;
  unique_ptr<double[]> coulomb(new double[max_coulomb_size]);
  map<pair<size_t, size_t>, shared_ptr<Matrix>> out;

  for (auto& spair : source_pvec) {
    // phase accumulated by moving an operator past the whole left ket block
    const int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%2) << 1);

    { // P_ab (x) I
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb-1);
      const BlockKey right_target = spair.right.key();

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // P_ab (x) I
        left_ops_->P_ab_copy_to(scratch.get(), spair.left.key(), i, j);
        kronecker_product_I_B(1.0, spair.right.nstates, false, tpair.left.nstates, spair.left.nstates, scratch.get(), tpair.left.nstates, out_block->data(), out_block->ndim());

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // - 1.0 <L'| A |L> (x) <R'| B |R>
      const BlockKey left_target(spair.left.nelea-1, spair.left.neleb);
      const BlockKey right_target(spair.right.nelea, spair.right.neleb-1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateAlpha}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateBeta}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = mo2e[(p+roffset)+(a+loffset)*norb+(j)*norb*norb+(i)*norb*norb*norb];
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, -1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // I (x) P_ab
      const BlockKey left_target = spair.left.key();
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb-1);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        // I (x) P_ab
        right_ops_->P_ab_copy_to(scratch.get(), spair.right.key(), i, j);
        kronecker_product_A_I(1.0, false, tpair.right.nstates, spair.right.nstates, scratch.get(), tpair.right.nstates, spair.left.nstates, out_block->data(), out_block->ndim());

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }

    { // + 1.0 <L'| B |L> (x) <R'| A |R>
      const BlockKey left_target(spair.left.nelea, spair.left.neleb-1);
      const BlockKey right_target(spair.right.nelea-1, spair.right.neleb);

      auto iter = find_if(target_pvec.begin(), target_pvec.end(), [&left_target, &right_target] (DMRG::BlockPair bp)
        { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }
      );
      if(iter!=target_pvec.end()) {
        DMRG::BlockPair tpair = *iter;

        auto out_block = make_shared<Matrix>(tpair.nstates(), spair.nstates(), true);

        shared_ptr<const btas::Tensor3<double>> Lgamma = blocks_->left_block()->coupling({GammaSQ::CreateBeta}).at({spair.left.key(),tpair.left.key()}).data;
        shared_ptr<const btas::Tensor3<double>> Rgamma = blocks_->right_block()->coupling({GammaSQ::CreateAlpha}).at({spair.right.key(),tpair.right.key()}).data;

        const int Lndim = Lgamma->extent(0);
        const int Lmdim = Lgamma->extent(1);
        const int Lsize = Lgamma->extent(0) * Lgamma->extent(1);

        const int Rndim = Rgamma->extent(0);
        const int Rmdim = Rgamma->extent(1);
        const int Rsize = Rgamma->extent(0) * Rgamma->extent(1);

        // fill in coulomb part
        assert(max_coulomb_size >= rnorb*lnorb);
        fill_n(scratch.get(), Lsize * rnorb, 0.0);
        for (int p = 0; p < rnorb; ++p) {
          for (int a = 0; a < lnorb; ++a) {
            coulomb[a + p*lnorb] = mo2e[(a+loffset)+(p+roffset)*norb+(j)*norb*norb+(i)*norb*norb*norb];
          }
        }

        dgemm_("N", "N", Lsize, rnorb, lnorb, 1.0, Lgamma->data(), Lsize, coulomb.get(), lnorb, 1.0, scratch.get(), Lsize);
        for (int p = 0; p < rnorb; ++p) {
          const double* ldata = scratch.get() + Lsize * p;
          const double* rdata = Rgamma->data() + Rsize * p;
          kronecker_product(left_phase, true, Rndim, Rmdim, rdata, Rndim, true, Lndim, Lmdim, ldata, Lndim, out_block->data(), out_block->ndim());
        }

        // add to map if large enough
        if (out_block->rms() > thresh_)
          out.emplace(make_pair<size_t, size_t>(tpair.offset, spair.offset), out_block);
      }
    }
  }

  return make_shared<BlockSparseMatrix>(target_info.nstates, source_info.nstates, out);
}

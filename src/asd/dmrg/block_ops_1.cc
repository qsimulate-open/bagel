//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg/block_ops_1.cc
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

using namespace std;
using namespace bagel;

BlockOperators1::BlockOperators1(shared_ptr<const DMRG_Block1> left, shared_ptr<DimerJop> jop, const double thresh) : BlockOperators(thresh), left_(left) {
  const int rnorb = jop->monomer_jop<0>()->nocc();
  const int lnorb = jop->monomer_jop<1>()->nocc();

  for (auto& binfo : left->blocks()) {
    const BlockKey bk = binfo.key();
    { // build pure parts
      shared_ptr<const btas::Tensor3<double>> gamma_aa = left->coupling({GammaSQ::CreateAlpha,GammaSQ::AnnihilateAlpha}).at({bk,bk}).data;
      shared_ptr<const btas::Tensor3<double>> gamma_bb = left->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta }).at({bk,bk}).data;

      const int gsize = gamma_aa->extent(0) * gamma_aa->extent(1);
      assert(gsize == binfo.nstates*binfo.nstates);

      // 2e part
      shared_ptr<Matrix> ham = left->h2e(bk)->copy();

      // 1e part
      const btas::Tensor3<double> gamma_aa_plus_bb(*gamma_aa + *gamma_bb);
      assert(ham->size()==gsize);
      shared_ptr<const Matrix> mo1e = jop->monomer_jop<1>()->mo1e()->matrix();
      dgemv_("N", gsize, mo1e->size(), 1.0, gamma_aa_plus_bb.data(), gsize, mo1e->data(), 1, 1.0, ham->data(), 1);
      ham_.emplace(bk, ham);

      // now for Q parts
      const Matrix& coulomb = *jop->coulomb_matrix<1,0,1,0>();
      const Matrix& exchange = *jop->coulomb_matrix<0,1,1,0>();

      auto Qaa = make_shared<btas::Tensor4<double>>(gamma_aa->extent(0), gamma_aa->extent(1), rnorb, rnorb);
      dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb,  1.0, gamma_aa_plus_bb.data(), gsize, coulomb.data(), coulomb.ndim(),
                                                         0.0, Qaa->data(), gsize);
      dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb, -1.0, gamma_aa->data(), gsize, exchange.data(), exchange.ndim(),
                                                         1.0, Qaa->data(), gsize);
      Q_aa_.emplace(bk, Qaa);

      auto Qbb = make_shared<btas::Tensor4<double>>(gamma_bb->extent(0), gamma_bb->extent(1), rnorb, rnorb);
      dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb,  1.0, gamma_aa_plus_bb.data(), gsize, coulomb.data(), coulomb.ndim(),
                                                         0.0, Qbb->data(), gsize);
      dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb, -1.0, gamma_bb->data(), gsize, exchange.data(), exchange.ndim(),
                                                         1.0, Qbb->data(), gsize);
      Q_bb_.emplace(bk, Qbb);

      const BlockKey abkey(bk.nelea-1, bk.neleb+1);
      if (left->contains(abkey)) {
        shared_ptr<btas::Tensor3<double>> gamma_ab = left->coupling({GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}).at({abkey,bk}).data;
        auto Qab = make_shared<btas::Tensor4<double>>(gamma_ab->extent(0), gamma_ab->extent(1), rnorb, rnorb);
        const int gabsize = gamma_ab->extent(0)*gamma_ab->extent(1);
        dgemm_("N", "T", gabsize, rnorb*rnorb, lnorb*lnorb, -1.0, gamma_ab->data(), gabsize, exchange.data(), exchange.ndim(),
                                                             0.0, Qab->data(), gabsize);
        Q_ab_.emplace(bk, Qab);
      }
    }

    { // S_a and S_b
      const Matrix& J_1101 = *jop->coulomb_matrix<1,1,0,1>();
      const Matrix& h01 = *jop->cross_mo1e();

      const BlockKey akey(bk.nelea+1, bk.neleb);
      if (left->contains(akey)) {
        shared_ptr<const btas::Tensor3<double>> gamma_a = left->coupling({GammaSQ::CreateAlpha}).at({akey,bk}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_aaa = left->coupling({GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({akey,bk}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_abb = left->coupling({GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({akey,bk}).data;

        const int asize = gamma_aaa->extent(0) * gamma_aaa->extent(1);

        // S_a
        btas::Tensor3<double> gamma_akk(*gamma_aaa + *gamma_abb);
        auto Sa = make_shared<btas::Tensor3<double>>(gamma_akk.extent(0), gamma_akk.extent(1), rnorb);
        dgemm_("N", "T", asize, rnorb,             lnorb, 1.0,  gamma_a->data(), asize,    h01.data(),    h01.ndim(), 0.0, Sa->data(), asize);
        dgemm_("N", "T", asize, rnorb, lnorb*lnorb*lnorb, 1.0, gamma_akk.data(), asize, J_1101.data(), J_1101.ndim(), 1.0, Sa->data(), asize);
        S_a_.emplace(bk, Sa);
      }

      const BlockKey bkey(bk.nelea, bk.neleb+1);
      if (left->contains(bkey)) {
        shared_ptr<const btas::Tensor3<double>> gamma_b = left->coupling({GammaSQ::CreateBeta}).at({bkey,bk}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_baa = left->coupling({GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}).at({bkey, bk}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_bbb = left->coupling({GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}).at({bkey, bk}).data;

        const int bsize = gamma_bbb->extent(0) * gamma_bbb->extent(1);

        // S_b
        btas::Tensor3<double> gamma_bkk(*gamma_baa + *gamma_bbb);
        auto Sb = make_shared<btas::Tensor3<double>>(gamma_bkk.extent(0), gamma_bkk.extent(1), rnorb);
        dgemm_("N", "T", bsize, rnorb,             lnorb, 1.0,  gamma_b->data(), bsize,    h01.data(),    h01.ndim(), 0.0, Sb->data(), bsize);
        dgemm_("N", "T", bsize, rnorb, lnorb*lnorb*lnorb, 1.0, gamma_bkk.data(), bsize, J_1101.data(), J_1101.ndim(), 1.0, Sb->data(), bsize);
        S_b_.emplace(bk, Sb);
      }
    }

    { // P_aa, P_bb, P_ab
      const Matrix& J_0110 = *jop->coulomb_matrix<0,1,1,0>();
      auto compute_Pxx = [&J_0110, &lnorb, &rnorb, &left, &bk] (const BlockKey new_key, list<GammaSQ> ops, const double fac) {
        shared_ptr<const btas::Tensor4<double>> gamma = left->coupling(ops).at({new_key, bk}).data;
        auto Pxx = make_shared<btas::Tensor4<double>>(gamma->extent(0), gamma->extent(1), rnorb, rnorb);
        const int gsize = gamma->extent(0)*gamma->extent(1);
        dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb, fac, gamma->data(), gsize, J_0110.data(), J_0110.ndim(), 0.0, Pxx->data(), gsize);
        return Pxx;
      };

      const BlockKey aakey(bk.nelea-2,bk.neleb);
      if (left->contains(aakey))
        P_aa_.emplace(bk, compute_Pxx(aakey, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, 0.5));

      const BlockKey bbkey(bk.nelea,bk.neleb-2);
      if (left->contains(bbkey))
        P_bb_.emplace(bk, compute_Pxx(bbkey, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}, 0.5));

      const BlockKey abkey(bk.nelea-1, bk.neleb-1);
      if (left->contains(abkey))
        P_ab_.emplace(bk, compute_Pxx(abkey, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha}, 1.0));
    }

    { // D_a, D_b
      const Matrix& J_0100 = *jop->coulomb_matrix<0,1,0,0>();

      const BlockKey akey(bk.nelea+1, bk.neleb);
      if (left->contains(akey)) {
        shared_ptr<btas::Tensor3<double>> gamma_a = left->coupling({GammaSQ::CreateAlpha}).at({akey,bk}).data;
        auto Da = make_shared<btas::TensorN<double,5>>(gamma_a->extent(0), gamma_a->extent(1), rnorb, rnorb, rnorb);
        const int gsize = Da->extent(0)*Da->extent(1);
        dgemm_("N", "T", gsize, rnorb*rnorb*rnorb, lnorb, 1.0, gamma_a->data(), gsize, J_0100.data(), J_0100.ndim(), 0.0, Da->data(), gsize);
        D_a_.emplace(bk, Da);
      }

      const BlockKey bkey(bk.nelea, bk.neleb+1);
      if (left->contains(bkey)) {
        shared_ptr<btas::Tensor3<double>> gamma_b = left->coupling({GammaSQ::CreateBeta}).at({bkey,bk}).data;
        auto Db = make_shared<btas::TensorN<double,5>>(gamma_b->extent(0), gamma_b->extent(1), rnorb, rnorb, rnorb);
        const int gsize = Db->extent(0)*Db->extent(1);
        dgemm_("N", "T", gsize, rnorb*rnorb*rnorb, lnorb, 1.0, gamma_b->data(), gsize, J_0100.data(), J_0100.ndim(), 0.0, Db->data(), gsize);
        D_b_.emplace(bk, Db);
      }
    }
  }
}

shared_ptr<BlockSparseMatrix> BlockOperators1::gamma_a(const BlockKey bk, int i) const {
  return get_sparse_mat_block(left_->coupling({GammaSQ::CreateAlpha}).at({BlockKey(bk.nelea+1,bk.neleb), bk}).data, i);
}

shared_ptr<BlockSparseMatrix> BlockOperators1::gamma_b(const BlockKey bk, int i) const {
  return get_sparse_mat_block(left_->coupling({GammaSQ::CreateBeta}).at({BlockKey(bk.nelea,bk.neleb+1), bk}).data, i);
}

shared_ptr<Matrix> BlockOperators1::gamma_a_as_matrix(const BlockKey bk, int i) const {
  return get_mat_block(left_->coupling({GammaSQ::CreateAlpha}).at({BlockKey(bk.nelea+1,bk.neleb), bk}).data, i);
}

shared_ptr<Matrix> BlockOperators1::gamma_b_as_matrix(const BlockKey bk, int i) const {
  return get_mat_block(left_->coupling({GammaSQ::CreateBeta}).at({BlockKey(bk.nelea,bk.neleb+1), bk}).data, i);
}

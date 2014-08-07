//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg/block_operators.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd_dmrg/block_operators.h>

using namespace std;
using namespace bagel;

BlockOperators::BlockOperators(shared_ptr<const DMRG_Block> left, shared_ptr<DimerJop> jop) {
  const int lnorb = jop->monomer_jop<0>()->nocc();
  const int rnorb = jop->monomer_jop<1>()->nocc();

  for (auto& binfo : left->blocks()) {
    const BlockKey bk = binfo.key();
    { // build pure parts
      shared_ptr<const btas::Tensor3<double>> gamma_aa = left->coupling({GammaSQ::AnnihilateAlpha,GammaSQ::CreateAlpha}).at({bk,bk}).data;
      shared_ptr<const btas::Tensor3<double>> gamma_bb = left->coupling({GammaSQ::AnnihilateBeta,GammaSQ::CreateBeta}).at({bk,bk}).data;

      const int gsize = gamma_aa->extent(0) * gamma_aa->extent(1);

      // 2e part
      shared_ptr<Matrix> ham = left->h2e(bk)->copy();

      // 1e part
      btas::Tensor3<double> gamma_aa_plus_bb(*gamma_aa + *gamma_bb);
      assert(ham->size()==gsize);
      shared_ptr<const Matrix> mo1e = jop->monomer_jop<1>()->mo1e()->matrix();
      dgemv_("N", gsize, mo1e->size(), 1.0, gamma_aa_plus_bb.data(), gsize, mo1e->data(), 1, 1.0, ham->data(), 1);
      ham_.emplace(bk, ham);

      // now for Q parts
      const Matrix& coulomb = *jop->coulomb_matrix<1,0,1,0>();
      const Matrix& exchange = *jop->coulomb_matrix<1,1,0,0>();

      auto Qaa = make_shared<btas::Tensor4<double>>(gamma_aa->extent(0), gamma_aa->extent(1), rnorb, rnorb);
      dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb,  1.0, gamma_aa_plus_bb.data(), gsize, coulomb.data(), coulomb.ndim(),
                                                         0.0, Qaa->data(), gsize);
      dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb, -1.0, gamma_bb->data(), gsize, exchange.data(), exchange.ndim(),
                                                         1.0, Qaa->data(), gsize);
      Q_aa_.emplace(bk, Qaa);

      auto Qbb = make_shared<btas::Tensor4<double>>(gamma_bb->extent(0), gamma_bb->extent(1), rnorb, rnorb);
      dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb,  1.0, gamma_aa_plus_bb.data(), gsize, coulomb.data(), coulomb.ndim(),
                                                         0.0, Qbb->data(), gsize);
      dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb, -1.0, gamma_aa->data(), gsize, exchange.data(), exchange.ndim(),
                                                         1.0, Qbb->data(), gsize);
      Q_bb_.emplace(bk, Qbb);

      const BlockKey abkey(bk.nelea-1, bk.neleb+1);
      if (left->contains(abkey)) {
        shared_ptr<btas::Tensor3<double>> gamma_ab = left->coupling({GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta}).at({bk,abkey}).data;
        auto Qab = make_shared<btas::Tensor4<double>>(gamma_ab->extent(0), gamma_ab->extent(1), rnorb, rnorb);
        const int gabsize = gamma_ab->extent(0)*gamma_ab->extent(1);
        dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb, -1.0, gamma_ab->data(), gabsize, exchange.data(), exchange.ndim(),
                                                           0.0, Qab->data(), gabsize);
        Q_ab_.emplace(bk, Qab);
      }
    }

    { // S_a and S_b
      const Matrix& J_0111 = *jop->coulomb_matrix<0,1,1,1>();
      const Matrix& h01 = *jop->cross_mo1e();

      const BlockKey akey(bk.nelea+1, bk.neleb);
      if (left->contains(akey)) {
        shared_ptr<const btas::Tensor3<double>> gamma_a = left->coupling({GammaSQ::CreateAlpha}).at({bk,akey}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_aaa = left->coupling({GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}).at({bk, akey}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_abb = left->coupling({GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha}).at({bk, akey}).data;

        const int asize = gamma_aaa->extent(0) * gamma_aaa->extent(1);

        // S_a
        btas::Tensor3<double> gamma_akk(*gamma_aaa + *gamma_abb);
        auto Sa = make_shared<btas::Tensor3<double>>(gamma_akk.extent(0), gamma_akk.extent(1), rnorb);
        dgemm_("N", "T", asize, rnorb,             lnorb, 1.0,  gamma_a->data(), asize,    h01.data(),    h01.ndim(), 0.0, Sa->data(), asize);
        dgemm_("N", "T", asize, rnorb, lnorb*lnorb*lnorb, 1.0, gamma_akk.data(), asize, J_0111.data(), J_0111.ndim(), 1.0, Sa->data(), asize);
        S_a_.emplace(bk, Sa);
      }

      const BlockKey bkey(bk.nelea, bk.neleb+1);
      if (left->contains(bkey)) {
        shared_ptr<const btas::Tensor3<double>> gamma_b = left->coupling({GammaSQ::CreateBeta}).at({bk,bkey}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_baa = left->coupling({GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta}).at({bk, bkey}).data;
        shared_ptr<const btas::Tensor3<double>> gamma_bbb = left->coupling({GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta}).at({bk, bkey}).data;

        const int bsize = gamma_bbb->extent(0) * gamma_bbb->extent(1);

        // S_b
        btas::Tensor3<double> gamma_bkk(*gamma_baa + *gamma_bbb);
        auto Sb = make_shared<btas::Tensor3<double>>(gamma_bkk.extent(0), gamma_bkk.extent(1), rnorb);
        dgemm_("N", "T", bsize, rnorb,             lnorb, 1.0,  gamma_b->data(), bsize,    h01.data(),    h01.ndim(), 0.0, Sb->data(), bsize);
        dgemm_("N", "T", bsize, rnorb, lnorb*lnorb*lnorb, 1.0, gamma_bkk.data(), bsize, J_0111.data(), J_0111.ndim(), 1.0, Sb->data(), bsize);
        S_b_.emplace(bk, Sb);
      }
    }

    { // P_aa, P_bb, P_ab
      const Matrix& J_0011 = *jop->coulomb_matrix<0,0,1,1>();
      auto compute_Pxx = [&J_0011, &lnorb, &rnorb, &left, &bk] (const BlockKey new_key, list<GammaSQ> ops) {
        shared_ptr<const btas::Tensor4<double>> gamma = left->coupling(ops).at({bk, new_key}).data;
        auto Pxx = make_shared<btas::Tensor4<double>>(gamma->extent(0), gamma->extent(1), rnorb, rnorb);
        const int gsize = gamma->extent(0)*gamma->extent(1);
        dgemm_("N", "T", gsize, rnorb*rnorb, lnorb*lnorb, 1.0, gamma->data(), gsize, J_0011.data(), J_0011.ndim(), 0.0, Pxx->data(), gsize);
        return Pxx;
      };

      const BlockKey aakey(bk.nelea-2,bk.neleb);
      if (left->contains(aakey))
        P_aa_.emplace(bk, compute_Pxx(aakey, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}));

      const BlockKey bbkey(bk.nelea,bk.neleb-2);
      if (left->contains(bbkey))
        P_bb_.emplace(bk, compute_Pxx(bbkey, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta}));

      const BlockKey abkey(bk.nelea-1, bk.neleb-1);
      if (left->contains(abkey))
        P_ab_.emplace(bk, compute_Pxx(abkey, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta}));
    }

    { // D_a, D_b
      const Matrix& J_0010 = *jop->coulomb_matrix<0,0,1,0>();

      const BlockKey akey(bk.nelea+1, bk.neleb);
      if (left->contains(akey)) {
        shared_ptr<btas::Tensor3<double>> gamma_a = left->coupling({GammaSQ::CreateAlpha}).at({bk, akey}).data;
        auto Da = make_shared<btas::TensorN<double,5>>(gamma_a->extent(0), gamma_a->extent(1), rnorb, rnorb, rnorb);
        const int gsize = Da->extent(0)*Da->extent(1);
        dgemm_("N", "T", gsize, rnorb*rnorb*rnorb, lnorb, 1.0, gamma_a->data(), gsize, J_0010.data(), J_0010.ndim(), 0.0, Da->data(), gsize);
        D_a_.emplace(bk, Da);
      }

      const BlockKey bkey(bk.nelea, bk.neleb+1);
      if (left->contains(bkey)) {
        shared_ptr<btas::Tensor3<double>> gamma_b = left->coupling({GammaSQ::CreateBeta}).at({bk, bkey}).data;
        auto Db = make_shared<btas::TensorN<double,5>>(gamma_b->extent(0), gamma_b->extent(1), rnorb, rnorb, rnorb);
        const int gsize = Db->extent(0)*Db->extent(1);
        dgemm_("N", "T", gsize, rnorb*rnorb*rnorb, lnorb, 1.0, gamma_b->data(), gsize, J_0010.data(), J_0010.ndim(), 0.0, Db->data(), gsize);
        D_b_.emplace(bk, Db);
      }
    }
  }
}


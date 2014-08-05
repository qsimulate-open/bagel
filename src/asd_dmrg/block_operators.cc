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
  // build pure parts
  for (auto& bk : left->blocks()) {

    // useful all over
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

    const int lnorb = jop->monomer_jop<0>()->nocc();
    const int rnorb = jop->monomer_jop<1>()->nocc();

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
  }
}


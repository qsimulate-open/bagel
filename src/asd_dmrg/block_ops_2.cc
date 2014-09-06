//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg/block_ops_2.cc
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
#include <src/asd_dmrg/dmrg_block.h>

using namespace std;
using namespace bagel;

BlockOperators2::BlockOperators2(shared_ptr<const DMRG_Block2> blocks, shared_ptr<DimerJop> jop) : blocks_(blocks), jop_(jop) {
  const int norb = jop->nocc();
  Matrix mo1e = *jop->mo1e()->matrix();
  btas::TensorView4<double> mo2e = btas::make_view(btas::CRange<4>(norb, norb, norb, norb), jop->mo2e()->storage());

  const int rasnorb = norb - (blocks->left_block()->norb() + blocks->right_block()->norb());
  const int lnorb = blocks->left_block()->norb();
  const int rnorb = blocks->right_block()->norb();

  { // make left jop to be used for left BlockOperators1
    const int leftsize = rasnorb+lnorb;
    auto left1e = make_shared<CSymMatrix>(mo1e.get_submatrix(0, 0, leftsize, leftsize));
    auto left2e = make_shared<Matrix>(leftsize*leftsize, leftsize*leftsize);
    auto low = {0, 0, 0, 0};
    auto up = {leftsize, leftsize, leftsize, leftsize};
    btas::TensorView4<double> mo2eslice = btas::make_view(mo2e.range().slice(low, up), mo2e.storage());
    copy(mo2eslice.begin(), mo2eslice.end(), left2e->begin());
    auto leftjop = make_shared<DimerJop>(rasnorb, lnorb, left1e, left2e);
    left_ops_ = make_shared<BlockOperators1>(blocks->left_block(), leftjop);
  }

  { // make right jop to be used for right BlockOperators1
    const int rightsize = rasnorb + rnorb;
    auto right1efull = mo1e.resize(rightsize, rightsize); // keeps the rasnorb x rasnorb terms in there
    right1efull->copy_block(0, rasnorb, rasnorb, rnorb, mo1e.get_submatrix(0, rasnorb+lnorb, rasnorb, rnorb));
    right1efull->copy_block(rasnorb, 0, rnorb, rasnorb, mo1e.get_submatrix(rasnorb+lnorb, 0, rnorb, rasnorb));
    right1efull->copy_block(rasnorb, rasnorb, rnorb, rnorb, mo1e.get_submatrix(rasnorb+lnorb, rasnorb+lnorb, rnorb, rnorb));
    auto right1e = make_shared<CSymMatrix>(right1efull);

    auto right2e = make_shared<Matrix>(rightsize*rightsize, rightsize*rightsize);
    btas::TensorView4<double> right2eView = btas::make_view(btas::CRange<4>(rightsize, rightsize, rightsize, rightsize), right2e->storage());
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

shared_ptr<Matrix> BlockOperators2::ham(const BlockKey bk) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::Q_aa(const BlockKey bk, int i, int j) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::Q_bb(const BlockKey bk, int i, int j) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::Q_ab(const BlockKey bk, int i, int j) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::S_a(const BlockKey bk, int i) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::S_b(const BlockKey bk, int i) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::P_aa(const BlockKey bk, int i, int j) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::P_bb(const BlockKey bk, int i, int j) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::P_ab(const BlockKey bk, int i, int j) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::gamma_a(const BlockKey bk, int i) const { return nullptr; }
shared_ptr<Matrix> BlockOperators2::gamma_b(const BlockKey bk, int i) const { return nullptr; }

double BlockOperators2::D_a(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const { return 0.0; }
double BlockOperators2::D_b(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const { return 0.0; }

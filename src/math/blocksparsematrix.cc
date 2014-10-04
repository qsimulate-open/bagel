//
// BAGEL - Parallel electron correlation program.
// Filename: blocksparsematrix.h
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <src/math/blocksparsematrix.h>

using namespace std;
using namespace bagel;

BlockSparseMatrix::BlockSparseMatrix(shared_ptr<Matrix> m) : ndim_(m->ndim()), mdim_(m->mdim()) {
  data_.emplace(make_pair(0,0), m);
}

BlockSparseMatrix::BlockSparseMatrix(const int n, const int m, map<pair<size_t, size_t>, shared_ptr<Matrix>> d) : data_(d), ndim_(n), mdim_(m) {
  assert(all_of(d.begin(), d.end(), [&n, &m] (pair<pair<size_t, size_t>, shared_ptr<Matrix>> p) {
    return (p.first.first + p.second->ndim() <= n && p.first.second + p.second->mdim() <= m);
  }));
}

// TODO: optimize by avoiding calls to element()
shared_ptr<VectorB> BlockSparseMatrix::diagonal() const {
  const size_t outdim = min(ndim_, mdim_);

  auto out = make_shared<VectorB>(outdim);
  for (size_t i = 0; i < outdim; ++i)
    (*out)(i) = element(i,i);

  return out;
}

void mat_block_multiply(const bool Atrans, const bool Btrans, const double alpha, const Matrix& A, const BlockSparseMatrix& B, const double beta, Matrix& C) {
  assert((Atrans ? A.mdim() : A.ndim())==C.ndim() &&
         (Btrans ? B.ndim() : B.mdim())==C.mdim() &&
         (Atrans ? A.ndim() : A.mdim())==(Btrans ? B.mdim() : B.ndim()));

  string At(Atrans ? "T" : "N");
  string Bt(Btrans ? "T" : "N");

  const int n = Atrans ? A.mdim() : A.ndim();

  for (auto block : B.data()) {
    const int bnstart = block.first.second;
    const int bmstart = block.first.second;
    shared_ptr<const Matrix> bmat = block.second;

    const int m = Btrans ? bmat->ndim() : bmat->mdim();
    const int k = Btrans ? bmat->mdim() : bmat->ndim();

    const int kstart = Btrans ? bmstart : bmstart;

    const double* adata = Atrans ? A.element_ptr(kstart, 0) : A.element_ptr(0, kstart);

    double* cdata = Btrans ? C.element_ptr(0, bnstart) : C.element_ptr(0, bmstart);

    dgemm_(At.c_str(), Bt.c_str(), n, m, k, alpha, adata, A.ndim(), bmat->data(), bmat->ndim(), beta, cdata, C.ndim());
  }
}

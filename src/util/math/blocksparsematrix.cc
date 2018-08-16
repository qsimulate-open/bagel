//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: blocksparsematrix.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/util/math/blocksparsematrix.h>

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

void bagel::mat_block_multiply(const bool Atrans, const bool Btrans, const double alpha, const Matrix& A, const BlockSparseMatrix& B, const double beta, Matrix& C) {
  assert((Atrans ? A.mdim() : A.ndim())==C.ndim() &&
         (Btrans ? B.ndim() : B.mdim())==C.mdim() &&
         (Atrans ? A.ndim() : A.mdim())==(Btrans ? B.mdim() : B.ndim()));

  string At(Atrans ? "T" : "N");
  string Bt(Btrans ? "T" : "N");

  const int n = Atrans ? A.mdim() : A.ndim();

  C.scale(beta);

  for (auto& block : B.data()) {
    const size_t bnstart = block.first.first;
    const size_t bmstart = block.first.second;
    shared_ptr<const Matrix> bmat = block.second;

    const int m = Btrans ? bmat->ndim() : bmat->mdim();
    const int k = Btrans ? bmat->mdim() : bmat->ndim();

    const int kstart = Btrans ? bmstart : bnstart;

    const double* adata = Atrans ? A.element_ptr(kstart, 0) : A.element_ptr(0, kstart);

    double* cdata = Btrans ? C.element_ptr(0, bnstart) : C.element_ptr(0, bmstart);

    dgemm_(At.c_str(), Bt.c_str(), n, m, k, alpha, adata, A.ndim(), bmat->data(), bmat->ndim(), 1.0, cdata, C.ndim());
  }
}

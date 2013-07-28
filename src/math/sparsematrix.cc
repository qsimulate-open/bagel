//
// BAGEL - Parallel electron correlation program.
// Filename: sparsematrix.cc
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <src/util/f77.h>
#include <src/math/matrix.h>
#include <src/math/sparsematrix.h>

using namespace bagel;
using namespace std;

//Constructors
SparseMatrix::SparseMatrix(const int n, const int m, const bool sym) : ndim_(n), mdim_(m), ndiag_(std::min(ndim_,mdim_)), symmetric_(sym) {
  assert( !(symmetric_ && (ndim_ != mdim_)) ); // if symmetric matrix, ndim_ == mdim_

  diagonal_ = unique_ptr<double[]>(new double[ndiag_]);
}

SparseMatrix::SparseMatrix(const SparseMatrix& o) : ndim_(o.ndim_), mdim_(o.mdim_), ndiag_(o.ndiag_), symmetric_(o.symmetric_) {
  diagonal_ = unique_ptr<double[]>(new double[ndiag_]);
  copy_n(o.diagonal(), ndiag_, diagonal_.get());

  for (auto& iele : offdiagonal_) { offdiagonal_.emplace_back(iele); }
}

//Scalar operations
void SparseMatrix::scale(const double& a) {
  dscal_(ndiag_, a, diagonal_.get(), 1);
  for (auto& iele : offdiagonal_) { iele.value *= a; }
}

SparseMatrix SparseMatrix::operator*(const double& a) const {
  SparseMatrix out(*this);
  out.scale(a);
  return out;
}

SparseMatrix SparseMatrix::operator/(const double& a) const {
  SparseMatrix out(*this);
  out.scale(1.0/a);
  return out;
}

// Matrix operations
Matrix SparseMatrix::multiply(const Matrix& o, const double eyediag) const {
  const int l = this->mdim();
  assert(l == o.ndim());
  const int m = o.mdim();
  const int n = this->ndim();

  Matrix out(n, m);

  for (int j = 0; j < m; ++j) {
    double* target = out.element_ptr(0,j);
    const double* source = o.element_ptr(0,j);
    for (int i = 0; i < n; ++i) {
      target[i] = (diagonal(i) + eyediag) * source[i];
    }
    // Avoiding putting an if statement inside a tight loop
    if (symmetric_) {
      for (auto& offdiag : offdiagonal_) {
        target[offdiag.i] += offdiag.value * source[offdiag.j];
        target[offdiag.j] += offdiag.value * source[offdiag.i];
      }
    }
    else {
      for (auto& offdiag : offdiagonal_) {
        target[offdiag.i] += offdiag.value * source[offdiag.j];
      }
    }
  }

  return out;
}

Matrix SparseMatrix::operator*(const Matrix& o) const {
  return multiply(o, 0.0);
}

Matrix SparseMatrix::operator+(const Matrix& o) const {
  assert( (this->ndim() == o.ndim()) && (this->mdim() == o.mdim()) );

  Matrix out(o);

  for(int i = 0; i < ndiag_; ++i) {
    out.element(i,i) += this->diagonal(i);
  }

  if (symmetric_) {
    for(auto& iele : offdiagonal_) {
      out.element(iele.i, iele.j) += iele.value;
      out.element(iele.j, iele.i) += iele.value;
    }
  }
  else {
    for(auto& iele : offdiagonal_) {
      out.element(iele.i, iele.j) += iele.value;
    }
  }

  return out;
}

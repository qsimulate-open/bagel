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

#include <src/math/sparsematrix.h>

using namespace bagel;
using namespace std;

// assumes data is given with 1-based indexing
SparseMatrix::SparseMatrix(const int n, const int m, const int size, double* data, int* cols, int* rind) :
  data_(unique_ptr<double[]>(new double[size])), cols_(unique_ptr<int[]>(new int[size])), rind_(unique_ptr<int[]>(new int[m+1])),
  ndim_(n), mdim_(m), size_(size)
{
  copy_n(data, size_, data_.get());
  copy_n(cols, size_, cols_.get());
  copy_n(rind, ndim_ + 1, rind_.get());
}

// assumes data is given with 1-based indexing
SparseMatrix::SparseMatrix(const int n, const int m, vector<double>& data, vector<int>& cols, vector<int>& rind) :
  SparseMatrix(n, m, data.size(), data.data(), cols.data(), rind.data())
{
  assert(data.size() == cols.size());
  assert(rind.size() == n + 1);
}

// assumes coords is given with 0-based indexing
SparseMatrix::SparseMatrix(const int n, const int m, map<pair<int, int>, double>& coords) : ndim_(n), mdim_(m), size_(coords.size())
{
  data_ = unique_ptr<double[]>(new double[size_]);
  cols_ = unique_ptr<int[]>(new int[size_]);
  rind_ = unique_ptr<int[]>(new int[ndim_ + 1]);

  int current_element = 0;
  int current_row = -1;

  for (auto& element : coords) {
    const int i = element.first.first;
    const int j = element.first.second;
    const double value = element.second;

    data_[current_element] = value;
    cols_[current_element] = j + 1; // for 1-based indexing
    
    if ( i >= current_row ) {
      while ( current_row < i )
        rind_[++current_row] = current_element + 1; // for 1-based indexing
    }
    else assert(false);

    ++current_element;
  }

  fill(rind_.get() + current_row + 1, rind_.get() + ndim_ + 1, size_ + 1);
}

SparseMatrix::SparseMatrix(const SparseMatrix& o) : SparseMatrix(o.ndim_, o.mdim_, o.size_, o.data_.get(), o.cols_.get(), o.rind_.get()) {}

SparseMatrix::SparseMatrix(SparseMatrix&& o) : data_(move(o.data_)), cols_(move(o.cols_)), rind_(move(o.rind_)), ndim_(o.ndim_), mdim_(o.mdim_), size_(o.size_) {}

//Scalar operations
void SparseMatrix::scale(const double& a) {
  transform(data_.get(), data_.get() + size_, data_.get(), [&a] (double p) { return p * a; });
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
Matrix SparseMatrix::operator*(const Matrix& o) const {
  const int l = this->mdim();
  assert(l == o.ndim());
  const int m = o.mdim();
  const int n = this->ndim();

  Matrix out(n, m);

  for (int j = 0; j < m; ++j) {
    double* target = out.element_ptr(0,j);
    const double* source = o.element_ptr(0,j);

    for (int i = 0; i < n; ++i) {
      for (int rowdata = rind_[i] - 1; rowdata < rind_[i+1] - 1; ++rowdata) {
        target[i] += data_[rowdata] * source[cols_[rowdata] - 1];
      }
    }
  }

  return out;
}

Matrix SparseMatrix::operator+(const Matrix& o) const {
  assert( (this->ndim() == o.ndim()) && (this->mdim() == o.mdim()) );

  Matrix out(o);

  for (int i = 0; i < ndim_; ++i) {
    for (int rowdata = rind_[i] - 1; rowdata < rind_[i+1] - 1; ++rowdata) {
      out.element(i, cols_[rowdata] - 1) += data_[rowdata];
    }
  }

  return out;
}

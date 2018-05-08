//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: sparsematrix.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <iostream>
#include <iomanip>

#include <src/util/math/sparsematrix.h>
#include <src/util/math/algo.h>

using namespace bagel;
using namespace std;

// assumes data is given with 1-based indexing
SparseMatrix::SparseMatrix(const int n, const int m, const int size, const double* data, int* cols, int* rind) :
  data_(unique_ptr<double[]>(new double[size])), cols_(unique_ptr<int[]>(new int[size])), rind_(unique_ptr<int[]>(new int[n+1])),
  ndim_(n), mdim_(m), size_(size)
{
  copy_n(data, size_, data_.get());
  copy_n(cols, size_, cols_.get());
  copy_n(rind, ndim_ + 1, rind_.get());
}

// assumes data is given with 1-based indexing
SparseMatrix::SparseMatrix(const int n, const int m, const vector<double>& data, vector<int>& cols, vector<int>& rind) :
  SparseMatrix(n, m, data.size(), data.data(), cols.data(), rind.data())
{
  assert(data.size() == cols.size());
  assert(rind.size() == n + 1);
}

// assumes coords is given with 0-based indexing
SparseMatrix::SparseMatrix(const int n, const int m, const map<pair<int, int>, double>& coords) : ndim_(n), mdim_(m), size_(coords.size())
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

// assumes coords is given with 0-based indexing
SparseMatrix::SparseMatrix(const int n, const int m, const vector<tuple<int, int, double>>& coords) : ndim_(n), mdim_(m), size_(coords.size())
{
  data_ = unique_ptr<double[]>(new double[size_]);
  cols_ = unique_ptr<int[]>(new int[size_]);
  rind_ = unique_ptr<int[]>(new int[ndim_ + 1]);

  int current_element = 0;
  int current_row = -1;

  for (auto& element : coords) {
    const int i = get<0>(element);
    const int j = get<1>(element);
    const double value = get<2>(element);

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
  for_each(data_.get(), data_.get() + size_, [&a] (double& p) { p *= a; });
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
  dcsrmm_("N", n, m, l, 1.0, data_.get(), cols_.get(), rind_.get(), o.data(), o.ndim(), 0.0, out.data(), out.ndim());

  return out;
}

Matrix SparseMatrix::operator+(const Matrix& o) const {
  assert( (this->ndim() == o.ndim()) && (this->mdim() == o.mdim()) );

  Matrix out(o);

  for (int i = 0; i < ndim_; ++i) {
    for (int rowdata = rind_[i] - 1; rowdata < rind_[i+1] - 1; ++rowdata) {
      out(i, cols_[rowdata] - 1) += data_[rowdata];
    }
  }

  return out;
}

void SparseMatrix::print_block_structure(const size_t nsize, const size_t msize) const {
  assert( nsize * msize > 0 );
  const size_t nblocks = (ndim_ - 1) / nsize + 1;
  const size_t mblocks = (mdim_ - 1) / msize + 1;

  vector<bool> structure(nblocks * mblocks, false);

  for (int i = 0; i < ndim_; ++i) {
    for (int rowdata = rind_[i] - 1; rowdata < rind_[i+1] - 1; ++rowdata) {
      const size_t j = cols_[rowdata] - 1;
      const size_t si = i / nsize;
      const size_t sj = j / msize;
      structure[si + nblocks * sj] = true;
    }
  }

  for (size_t is = 0; is < nblocks; ++is) {
    for (size_t js = 0; js < mblocks; ++js) {
      cout << ( structure[is + nblocks * js] ? "  1" : "  0" );
    }
    cout << endl;
  }
}

void SparseMatrix::print_table() const {
  for (int i = 0; i < ndim_; ++i) {
    for (int rowdata = rind_[i] - 1; rowdata < rind_[i+1] - 1; ++rowdata) {
      cout << setw(6) << setprecision(0) << "(" <<  i << ", " << cols_[rowdata] - 1 << ")" << setw(20) << setprecision(8) << data_[rowdata] << endl;
    }
  }
}

shared_ptr<Matrix> SparseMatrix::matrix() const {
  auto out = make_shared<Matrix>(ndim_, mdim_);
  for (int i = 0; i < ndim_; ++i) {
    for (int rowdata = rind_[i] - 1; rowdata < rind_[i+1] - 1; ++rowdata)
      out->element(i, cols_[rowdata]-1) = data_[rowdata];
  }

  return out;
}

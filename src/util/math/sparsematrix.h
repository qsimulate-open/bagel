//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: sparsematrix.h
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

#ifndef __math_sparse_matrix_h
#define __math_sparse_matrix_h

#include <src/util/math/matrix.h>

// Designed to be compliant with MKL_SPARSE_BLAS in which CSR seems to be the most commonly used format
//  see the MKL documentation for a more detailed description of the storage formats and some of the capabilities

// Apparently, when MKL Sparse includes a dense matrix, it uses the indexing scheme to determine whether the dense
// matrix is row-major of column-major: row-major for 0-based indexing and column-major for 1-based indexing.
// Our matrices are column-major, thus we need to use 1-based indexing.

namespace bagel {

class SparseMatrix {
  protected:
    // These are standard for CSR storage
    std::unique_ptr<double[]> data_; // contiguous list of non-zero elements: length = size_
    std::unique_ptr<int[]>    cols_; // contiguous list of column id of each element: length = size_
    std::unique_ptr<int[]>    rind_; // the elements of row i can be found starting at data_[rind_[i]]
                                        // and ending (without including) at data_[rind_[i+1]]

    const int ndim_;
    const int mdim_;
    const int size_;

  public:
    SparseMatrix(const int n, const int m, const int size, const double* data, int* cols, int* rind);
    SparseMatrix(const int n, const int m, const std::vector<double>& data, std::vector<int>& cols, std::vector<int>& rind);
    SparseMatrix(const int n, const int m, const std::map<std::pair<int, int>, double>& coords);
    SparseMatrix(const int n, const int m, const std::vector<std::tuple<int, int, double>>& coords);
    SparseMatrix(const SparseMatrix& o); // copy constructor
    SparseMatrix(SparseMatrix&& o);      // move constructor

    // getting info
    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }

    int* cols() { return cols_.get(); }
    int* rind() { return rind_.get(); }

    int ndim() const { return ndim_; }
    int mdim() const { return mdim_; }
    int size() const { return size_; }

    // Scalar operations (only bare minimum right now)
    void scale(const double& a);
    SparseMatrix operator*(const double& a) const;
    SparseMatrix operator/(const double& a) const;

    // Matrix operations (bare minimum)
    Matrix multiply(const Matrix& o, const double eyediag) const;
    Matrix operator*(const Matrix& o) const;
    Matrix operator+(const Matrix& o) const;

    void zero() { std::fill_n(data_.get(), size_, 0.0); }

    // diagnostics
    void print_block_structure(const size_t nsize, const size_t msize) const;
    void print_table() const;

    std::shared_ptr<Matrix> matrix() const;
};

}

#endif

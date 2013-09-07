//
// BAGEL - Parallel electron correlation program.
// Filename: sparsematrix.h
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

#ifndef __math_sparse_matrix_h
#define __math_sparse_matrix_h

#include <src/math/matrix.h>

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
    SparseMatrix(const int n, const int m, const int size, double* data, int* cols, int* rind);
    SparseMatrix(const int n, const int m, std::vector<double>& data, std::vector<int>& cols, std::vector<int>& rind);
    SparseMatrix(const int n, const int m, std::map<std::pair<int, int>, double>& coords);
    SparseMatrix(const SparseMatrix& o); // copy constructor
    SparseMatrix(SparseMatrix&& o);      // move constructor

    // getting info
    const int ndim() const { return ndim_; }
    const int mdim() const { return mdim_; }
    const int size() const { return size_; }

    // Scalar operations (only bare minimum right now)
    void scale(const double& a);
    SparseMatrix operator*(const double& a) const;
    SparseMatrix operator/(const double& a) const;

    // Matrix operations (bare minimum)
    Matrix multiply(const Matrix& o, const double eyediag) const;
    Matrix operator*(const Matrix& o) const;
    Matrix operator+(const Matrix& o) const;
};

}

#endif

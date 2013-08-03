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

namespace bagel {

class SparseMatrix {
  protected:
    struct Element {
      const int i;
      const int j;
      double value;

      Element(const int ii, const int jj, const double vv) : i(ii), j(jj), value(vv) {}
      Element(const Element& o) : i(o.i), j(o.j), value(o.value) {}
    };

    const int ndim_;
    const int mdim_;
    const int ndiag_;

    const bool symmetric_;

    std::unique_ptr<double[]> diagonal_;
    std::vector<Element> offdiagonal_;

  public:
    // SparseMatrix needs to be built element-by-element, for the most part
    SparseMatrix(const int n, const int m, const bool sym = true);
    SparseMatrix(const SparseMatrix& o);

    // getting info
    const int ndim() const { return ndim_; }
    const int mdim() const { return mdim_; }
    const int ndiag() const { return ndiag_; }

    // Scalar operations (only bare minimum right now)
    void scale(const double& a);
    SparseMatrix operator*(const double& a) const;
    SparseMatrix operator/(const double& a) const;

    // Matrix operations (bare minimum)
    Matrix multiply(const Matrix& o, const double eyediag) const;
    Matrix operator*(const Matrix& o) const;
    Matrix operator+(const Matrix& o) const;

    double* diagonal() { return diagonal_.get(); }
    const double* diagonal() const { return diagonal_.get(); }
    double& diagonal(const int i) { return diagonal_[i]; }
    const double diagonal(const int i) const { return diagonal_[i]; }

    const std::vector<Element>& offdiagonal() const { return offdiagonal_; }

    void insert(const int i, const int j, const double val) { offdiagonal_.emplace_back(i,j,val); }
};

}

#endif

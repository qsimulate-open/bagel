//
// BAGEL - Parallel electron correlation program.
// Filename: matrix.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_UTIL_MATRIX_H
#define __SRC_UTIL_MATRIX_H

#include <string>
#include <list>
#include <src/util/f77.h>
#include <src/math/matrix_base.h>
#include <src/math/distmatrix_base.h>

namespace bagel {

#ifdef HAVE_SCALAPACK
class DistMatrix;
#else
class Matrix;
using DistMatrix = Matrix;
#endif

class Matrix : public Matrix_base<double>, public std::enable_shared_from_this<Matrix> {
  public:
#ifdef HAVE_SCALAPACK
    Matrix(const int n, const int m, const bool localized = false);
#else
    Matrix(const int n, const int m, const bool localized = true);
#endif
    Matrix(const Matrix&);

    std::shared_ptr<Matrix> cut(const int nstart, const int nend) const { return get_submatrix(nstart, 0, nend-nstart, mdim_); }
    std::shared_ptr<Matrix> slice(const int mstart, const int mend) const { return get_submatrix(0, mstart, ndim_, mend-mstart); }
    std::shared_ptr<Matrix> resize(const int n, const int m) const { return this->resize_impl<Matrix>(n, m); }
    std::shared_ptr<Matrix> merge(const std::shared_ptr<const Matrix> o) const { return this->merge_impl<Matrix>(o); }

    // antisymmetrize
    void antisymmetrize();

    // diagonalize this matrix (overwritten by a coefficient matrix)
    void diagonalize(double* vec) override;
    std::shared_ptr<Matrix> diagonalize_blocks(double* eig, std::vector<int> blocks) { return diagonalize_blocks_impl<Matrix>(eig, blocks); }
    void svd(std::shared_ptr<Matrix>, std::shared_ptr<Matrix>);
    // compute S^-1. Assumes positive definite matrix
    void inverse();
    // compute S^-1 using symmetric form.
    void inverse_symmetric(const double thresh = 1.0e-8);
    // compute S^-1/2. If an eigenvalue of S is smaller than thresh, the root will be discarded.
    void inverse_half(const double thresh = 1.0e-8);
    // compute S^1/2. Same algorithm as above.
    void sqrt();

    using Matrix_base<double>::copy_block;
    using Matrix_base<double>::get_block;
    using Matrix_base<double>::add_block;

    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, const Matrix& o) {
      assert(o.ndim() == nsize && o.mdim() == msize);
      this->copy_block(nstart, mstart, nsize, msize, o.data());
    }
    void add_block(const double a, const int nstart, const int mstart, const int nsize, const int msize, const Matrix& o) {
      assert(o.ndim() == nsize && o.mdim() == msize);
      this->add_block(a, nstart, mstart, nsize, msize, o.data());
    }

    std::shared_ptr<Matrix> get_submatrix(const int nstart, const int mstart, const int ndim, const int mdim) const {
      return this->get_submatrix_impl<Matrix>(nstart, mstart, ndim, mdim);
    }

    Matrix operator*(const Matrix&) const;
    Matrix& operator*=(const Matrix&);
    Matrix operator*(const double& a) const;
    Matrix operator/(const double& a) const;
    Matrix& operator*=(const double& a);
    Matrix& operator/=(const double& a);
    Matrix operator%(const Matrix&) const; // caution
    Matrix operator^(const Matrix&) const; // caution
    Matrix operator+(const Matrix&) const;
    Matrix& operator+=(const Matrix&);
    Matrix& operator-=(const Matrix&);
    Matrix& operator=(const Matrix&);
    Matrix operator-(const Matrix&) const;

    Matrix& operator/=(const Matrix&);
    Matrix operator/(const Matrix&) const;

    std::shared_ptr<Matrix> clone() const { return std::make_shared<Matrix>(ndim_, mdim_, localized_); }
    std::shared_ptr<Matrix> copy() const { return std::make_shared<Matrix>(*this); }

    // returns exp(*this)
    std::shared_ptr<Matrix> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<Matrix> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<Matrix> transpose(const double a = 1.0) const;

    using Matrix_base<double>::ax_plus_y;
    using Matrix_base<double>::dot_product;
    void ax_plus_y(const double a, const Matrix& o) { this->ax_plus_y_impl(a, o); }
    double dot_product(const Matrix& o) const { return this->dot_product_impl(o); }

    double orthog(const std::list<std::shared_ptr<const Matrix>> o) { return this->orthog_impl(o); } 
    void rotate(const int i, const int j, const double c, const double s) { drot_(ndim_, element_ptr(0,i), 1, element_ptr(0,j), 1, c, s); }
    void rotate(const int i, const int j, const double gamma) { rotate(i, j, cos(gamma), sin(gamma)); }

    // purify a (near unitary) matrix to be unitary
    void purify_unitary();
    void purify_idempotent(const Matrix& s);
    void purify_redrotation(const int nclosed, const int nact, const int nvirt);

    std::shared_ptr<Matrix> solve(std::shared_ptr<const Matrix> A, const int n) const;

    virtual void print(const std::string in = "", const size_t size = 10) const;

#ifdef HAVE_SCALAPACK
    // return a shared pointer to this ifndef HAVE_SCALAPACK
    std::shared_ptr<DistMatrix> distmatrix() const;

    Matrix(const DistMatrix&);
#else
    std::shared_ptr<const Matrix> distmatrix() const;
    std::shared_ptr<const Matrix> matrix() const { return shared_from_this(); }
    std::shared_ptr<const Matrix> form_density_rhf(const int n, const int off = 0) const;
#endif
};


#ifdef HAVE_SCALAPACK
// Not to be confused with Matrix. DistMatrix is distributed and only supported when SCALAPACK is turned on. Limited functionality
class DistMatrix : public DistMatrix_base<double> {
  public:
    DistMatrix(const int n, const int m);
    DistMatrix(const DistMatrix&);
    DistMatrix(const Matrix&);

    void diagonalize(double* vec) override;

    DistMatrix operator*(const DistMatrix&) const;
    DistMatrix& operator*=(const DistMatrix&);
    DistMatrix operator%(const DistMatrix&) const; // caution
    DistMatrix operator^(const DistMatrix&) const; // caution
    DistMatrix operator+(const DistMatrix& o) const { DistMatrix out(*this); out.ax_plus_y(1.0, o); return out; }
    DistMatrix operator-(const DistMatrix& o) const { DistMatrix out(*this); out.ax_plus_y(-1.0, o); return out; }
    DistMatrix& operator+=(const DistMatrix& o) { ax_plus_y(1.0, o); return *this; }
    DistMatrix& operator-=(const DistMatrix& o) { ax_plus_y(-1.0, o); return *this; }
    DistMatrix& operator=(const DistMatrix& o) { assert(size() == o.size()); std::copy_n(o.local_.get(), size(), local_.get()); return *this; }

    std::shared_ptr<DistMatrix> clone() const { return std::make_shared<DistMatrix>(ndim_, mdim_); }

    using DistMatrix_base<double>::scale;
    using DistMatrix_base<double>::ax_plus_y;
    using DistMatrix_base<double>::dot_product;

    void ax_plus_y(const double a, const DistMatrix& o) { this->ax_plus_y_impl(a,o); }
    double dot_product(const DistMatrix& o) const { return this->dot_product_impl(o); }

    void rotate(const int i, const int j, const double cs, const double sn);
    void rotate(const int i, const int j, const double gamma) { rotate(i, j, cos(gamma), sin(gamma)); }

    std::shared_ptr<Matrix> matrix() const;

    std::shared_ptr<const DistMatrix> form_density_rhf(const int n, const int off = 0) const;
};
#endif

}

#endif

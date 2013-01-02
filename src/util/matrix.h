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


#ifndef __SRC_UTIL_MATRIX_H
#define __SRC_UTIL_MATRIX_H

#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <config.h>
#include <src/util/f77.h>
#include <src/util/matrix_base.h>
#include <src/util/distmatrix_base.h>

namespace bagel {

#ifdef HAVE_SCALAPACK
class DistMatrix;
#else
class Matrix;
typedef Matrix DistMatrix;
#endif

class Matrix : public Matrix_base<double>, public std::enable_shared_from_this<Matrix> {
  public:
    Matrix(const int n, const int m);
    Matrix(const Matrix&);

    void antisymmetrize();
    std::shared_ptr<Matrix> cut(const int) const;
    std::shared_ptr<Matrix> resize(const int, const int) const;
    std::shared_ptr<Matrix> slice(const int, const int) const;
    std::shared_ptr<Matrix> merge(const std::shared_ptr<const Matrix>) const;
    // diagonalize this matrix (overwritten by a coefficient matrix)
    void diagonalize(double* vec) override;
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

    void copy_block(const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const Matrix> o);
    std::shared_ptr<Matrix> get_submatrix(const int nstart, const int mstart, const int ndim, const int mdim) const;
    void add_block(const int nstart, const int mstart, const int ndim, const int mdim, const Matrix& o);
    void add_block(const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const Matrix> o) { add_block(nstart, mstart, ndim, mdim, *o); };
    void add_block(const int nstart, const int mstart, const int nsize, const int msize, const double* o);

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

    std::shared_ptr<Matrix> clone() const { return std::shared_ptr<Matrix>(new Matrix(ndim_, mdim_)); }

    // returns exp(*this)
    std::shared_ptr<Matrix> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<Matrix> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<Matrix> transpose() const;

    void daxpy(const double, const Matrix&);
    void daxpy(const double, const std::shared_ptr<const Matrix>);
    double ddot(const Matrix&) const;
    double norm() const { return std::sqrt(ddot(*this)); }
    double ddot(const std::shared_ptr<const Matrix>) const;
    double rms() const;
    double trace() const;

    void dscal(const double a) { dscal_(size(), a, data(), 1); }
    void scale(const double a) { dscal(a); }

    void add_diag(const double a, const int i, const int j) {
      assert(ndim_ == mdim_);
      for (int ii = i; ii != j; ++ii) element(ii,ii) += a;
    }
    void add_diag(const double a) { add_diag(a,0,ndim_); }
    // returns diagonal elements
    std::unique_ptr<double[]> diag() const;

    void unit() { fill(0.0); for (int i = 0; i != ndim_; ++i) element(i,i) = 1.0; assert(ndim_ == mdim_);}
    // purify a (near unitary) matrix to be unitary

    void purify_unitary();
    void purify_idempotent(const Matrix& s);
    void purify_redrotation(const int nclosed, const int nact, const int nvirt);

    double orthog(const std::list<std::shared_ptr<const Matrix> > o);

    void print(const std::string in = "", const int size = 10) const;

    // return a shared pointer to this ifndef HAVE_SCALAPACK
    std::shared_ptr<const DistMatrix> distmatrix() const;

#ifdef HAVE_SCALAPACK
    Matrix(const DistMatrix&);
#else
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
    DistMatrix operator+(const DistMatrix& o) const { DistMatrix out(*this); out.daxpy(1.0, o); return out; }
    DistMatrix operator-(const DistMatrix& o) const { DistMatrix out(*this); out.daxpy(-1.0, o); return out; }
    DistMatrix& operator+=(const DistMatrix& o) { daxpy(1.0, o); return *this; }
    DistMatrix& operator-=(const DistMatrix& o) { daxpy(-1.0, o); return *this; }
    DistMatrix& operator=(const DistMatrix& o) { assert(size() == o.size()); std::copy_n(o.local_.get(), size(), local_.get()); return *this; }

    std::shared_ptr<DistMatrix> clone() const { return std::shared_ptr<DistMatrix>(new DistMatrix(ndim_, mdim_)); }

    void daxpy(const double a, const DistMatrix& o) { assert(size() == o.size()); daxpy_(size(), a, o.local_.get(), 1, local_.get(), 1); }
    void daxpy(const double a, const std::shared_ptr<const DistMatrix> o) { daxpy(a, *o); }

    double ddot(const DistMatrix&) const;
    double norm() const { return std::sqrt(ddot(*this)); }
    double ddot(const std::shared_ptr<const DistMatrix> o) const { return ddot(*o); }
    double rms() const { return norm()/std::sqrt(ndim_*mdim_); }

    void scale(const double a) { dscal_(size(), a, local_.get(), 1); }

    std::shared_ptr<Matrix> matrix() const;

    std::shared_ptr<const DistMatrix> form_density_rhf(const int n, const int off = 0) const;
};
#endif

}

#endif

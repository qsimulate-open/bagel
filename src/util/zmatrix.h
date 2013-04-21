//
// BAGEL - Parallel electron correlation program.
// Filename: zmatrix.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@u.northwestern.edu>
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


#ifndef __SRC_UTIL_ZMATRIX_H
#define __SRC_UTIL_ZMATRIX_H

#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <complex>
#include <src/util/matrix.h>

namespace bagel {

#ifdef HAVE_SCALAPACK
class DistZMatrix;
#else
class ZMatrix;
using DistZMatrix = ZMatrix;
#endif

class ZMatrix : public Matrix_base<std::complex<double>>, public std::enable_shared_from_this<ZMatrix> {
  public:
#ifdef HAVE_SCALAPACK
    ZMatrix(const int n, const int m, const bool localized = false);
#else
    ZMatrix(const int n, const int m, const bool localized = true);
#endif
    ZMatrix(const ZMatrix&);
    ZMatrix(const Matrix& real, const Matrix& imag);
    ZMatrix(const Matrix& real, const std::complex<double> factor);

    void antisymmetrize();
    void hermite();
    std::shared_ptr<ZMatrix> cut(const int) const;
    std::shared_ptr<ZMatrix> resize(const int, const int) const;
    std::shared_ptr<ZMatrix> slice(const int, const int) const;
    std::shared_ptr<ZMatrix> merge(const std::shared_ptr<const ZMatrix>) const;
    // diagonalize this matrix (overwritten by a coefficient matrix)
    virtual void diagonalize(double* vec);
    void diagonalize_skew(double* vec);

    void svd(std::shared_ptr<ZMatrix>, std::shared_ptr<ZMatrix>);
    // compute S^-1. Assumes positive definite matrix
    void inverse();
    // compute S^-1/2. If an eigenvalue of S is smaller than thresh, the root will be discarded.
    void inverse_half(const double thresh = 1.0e-8);

    using Matrix_base<std::complex<double>>::copy_block;
    using Matrix_base<std::complex<double>>::get_block;
    void copy_block(const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const ZMatrix> o);
    // copy a block from a Matrix object to ZMatrix object and mutliply by coefficient coeff. If type = 0, coeff is purely real, else purely imaginary
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const double* data);
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::unique_ptr<double[]> o);
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const Matrix> o);

    void add_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const double* data);
    void add_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::unique_ptr<double[]> o);
    void add_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const Matrix> o);

    void add_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::complex<double>* data);
    void add_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::unique_ptr<std::complex<double>[]> o);
    void add_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const ZMatrix> o);

    std::shared_ptr<Matrix> get_real_part() const;
    std::shared_ptr<Matrix> get_imag_part() const;

    std::shared_ptr<ZMatrix> get_conjg() const;

    std::shared_ptr<ZMatrix> get_submatrix(const int, const int, const int, const int) const;

    ZMatrix operator*(const ZMatrix&) const;
    ZMatrix& operator*=(const ZMatrix&);
    ZMatrix operator*(const std::complex<double>& a) const;
    ZMatrix operator/(const std::complex<double>& a) const;
    ZMatrix& operator*=(const std::complex<double>& a);
    ZMatrix& operator/=(const std::complex<double>& a);
    ZMatrix operator%(const ZMatrix&) const; // caution
    ZMatrix operator^(const ZMatrix&) const; // caution
    ZMatrix operator+(const ZMatrix&) const;
    ZMatrix& operator+=(const ZMatrix&);
    ZMatrix& operator-=(const ZMatrix&);
    ZMatrix& operator=(const ZMatrix&);
    ZMatrix operator-(const ZMatrix&) const;

    ZMatrix& operator/=(const ZMatrix&);
    ZMatrix operator/(const ZMatrix&) const;

    std::shared_ptr<ZMatrix> clone() const { return std::shared_ptr<ZMatrix>(new ZMatrix(ndim_, mdim_)); }
    std::shared_ptr<ZMatrix> copy() const { return std::shared_ptr<ZMatrix>(new ZMatrix(*this)); }

    // returns exp(*this)
    std::shared_ptr<ZMatrix> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<ZMatrix> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<ZMatrix> transpose() const;
    // returns hermite-conjugate(*this)
    std::shared_ptr<ZMatrix> transpose_conjg() const;

    void zaxpy(const std::complex<double>, const ZMatrix&);
    void zaxpy(const std::complex<double>, const std::shared_ptr<const ZMatrix>);
    std::complex<double> zdotc(const ZMatrix&) const;
    std::complex<double> zdotc(const std::shared_ptr<const ZMatrix>) const;
    double norm() const;
    double rms() const;
    std::complex<double> trace() const;

    // for generic algorithms 
    void daxpy(const std::complex<double> a, const ZMatrix& o) { zaxpy(a, o); }
    void daxpy(const std::complex<double> a, const std::shared_ptr<const ZMatrix> o) { zaxpy(a, o); }

    void zscal(const std::complex<double> a) { zscal_(size(), a, data(), 1); }
    void scale(const std::complex<double> a) { zscal(a); }
    std::complex<double> ddot(const ZMatrix& o) const { return zdotc(o); }
    std::complex<double> ddot(const std::shared_ptr<const ZMatrix> o) const { return zdotc(o); }

    void add_diag(const std::complex<double> a, const int i, const int j) {
      assert(ndim_ == mdim_);
      for (int ii = i; ii != j; ++ii) element(ii,ii) += a;
    }
    void add_diag(const std::complex<double> a) { add_diag(a,0,ndim_); }
    // returns diagonal elements
    std::unique_ptr<std::complex<double>[]> diag() const;

    void unit() { fill(std::complex<double>(0.0, 0.0)); for (int i = 0; i != ndim_; ++i) element(i,i) = std::complex<double>(1.0, 0.0); assert(ndim_ == mdim_);}
    // purify a (near unitary) matrix to be unitary

    void purify_unitary();
    void purify_idempotent(const ZMatrix& s);
    void purify_redrotation(const int nclosed, const int nact, const int nvirt);

    std::complex<double> orthog(const std::list<std::shared_ptr<const ZMatrix>> o);

    std::shared_ptr<ZMatrix> solve(std::shared_ptr<const ZMatrix> A, const int n) const;

    // first parameter must be "R", "I", or "T" to print real part, imaginary part, or total
    void print(const std::string, const std::string in = "", const size_t size = 10) const;
    void print_row(const std::string, const std::string in = "", const size_t size = 10, const int row_num = 0) const;
    void print_col(const std::string, const std::string in = "", const size_t size = 10, const int col_num = 0) const;

    std::shared_ptr<ZMatrix> convert_real(const std::shared_ptr<const Matrix>) const;

#ifdef HAVE_SCALAPACK
    std::shared_ptr<DistZMatrix> distmatrix() const;
    ZMatrix(const DistZMatrix&);
#else
    std::shared_ptr<const ZMatrix> distmatrix() const;
    std::shared_ptr<const ZMatrix> matrix() const { return shared_from_this(); }
    std::shared_ptr<const ZMatrix> form_density_rhf(const int n, const int off = 0) const;
#endif
};


#ifdef HAVE_SCALAPACK
// Not to be confused with Matrix. DistMatrix is distributed and only supported when SCALAPACK is turned on. Limited functionality 
class DistZMatrix : public DistMatrix_base<std::complex<double>> {
  public:
    DistZMatrix(const int n, const int m);
    DistZMatrix(const DistZMatrix&);
    DistZMatrix(const ZMatrix&);

    void diagonalize(double* vec) override;

    DistZMatrix operator*(const DistZMatrix&) const;
    DistZMatrix& operator*=(const DistZMatrix&);
    DistZMatrix operator%(const DistZMatrix&) const; // caution
    DistZMatrix operator^(const DistZMatrix&) const; // caution
    DistZMatrix operator+(const DistZMatrix& o) const { DistZMatrix out(*this); out.zaxpy(1.0, o); return out; }
    DistZMatrix operator-(const DistZMatrix& o) const { DistZMatrix out(*this); out.zaxpy(-1.0, o); return out; }
    DistZMatrix& operator+=(const DistZMatrix& o) { zaxpy(1.0, o); return *this; }
    DistZMatrix& operator-=(const DistZMatrix& o) { zaxpy(-1.0, o); return *this; }
    DistZMatrix& operator=(const DistZMatrix& o) { assert(size() == o.size()); std::copy_n(o.local_.get(), size(), local_.get()); return *this; }

    std::shared_ptr<DistZMatrix> clone() const { return std::shared_ptr<DistZMatrix>(new DistZMatrix(ndim_, mdim_)); }

    void zaxpy(const double a, const DistZMatrix& o) { const std::complex<double> b(a); zaxpy(b,o); }
    void zaxpy(const double a, const std::shared_ptr<const DistZMatrix> o) { zaxpy(a,*o); }
    void zaxpy(const std::complex<double> a, const DistZMatrix& o) { assert(size() == o.size()); zaxpy_(size(), a, o.local_.get(), 1, local_.get(), 1); }
    void zaxpy(const std::complex<double> a, const std::shared_ptr<const DistZMatrix> o) { zaxpy(a, *o); }

    std::complex<double> zdotc(const DistZMatrix&) const;
    std::complex<double> zdotc(const std::shared_ptr<const DistZMatrix> o) const { return zdotc(*o); }
    double norm() const { return std::sqrt(zdotc(*this).real()); }
    double rms() const { return norm()/std::sqrt(ndim_*mdim_); }

    // for generic algorithms
    std::complex<double> ddot(const DistZMatrix& o) const { return zdotc(o); }
    std::complex<double> ddot(const std::shared_ptr<const DistZMatrix> o) const { return zdotc(o); }
    void daxpy(const std::complex<double> a, const DistZMatrix& o) { zaxpy(a, o); } 
    void daxpy(const std::complex<double> a, const std::shared_ptr<const DistZMatrix> o) { zaxpy(a, o); }

    void scale(const std::complex<double> a) { zscal_(size(), a, local_.get(), 1); }
    void scale(const double a) { const std::complex<double> b(a); scale(b); }

    std::shared_ptr<ZMatrix> matrix() const;

    std::shared_ptr<const DistZMatrix> form_density_rhf(const int n, const int off = 0) const;
};
#endif

}

#endif

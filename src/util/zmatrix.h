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
#include <src/util/f77.h>
#include <complex>
#include <src/util/matrix.h>

namespace bagel {

class ZMatrix : public Matrix_base<std::complex<double> > {
  public:
    ZMatrix(const int n, const int m);
    ZMatrix(const ZMatrix&);

    void antisymmetrize();
    std::shared_ptr<ZMatrix> cut(const int) const;
    std::shared_ptr<ZMatrix> resize(const int, const int) const;
    std::shared_ptr<ZMatrix> slice(const int, const int) const;
    std::shared_ptr<ZMatrix> merge(const std::shared_ptr<const ZMatrix>) const;
    // diagonalize this matrix (overwritten by a coefficient matrix)
    virtual void diagonalize(double* vec);
    void svd(std::shared_ptr<ZMatrix>, std::shared_ptr<ZMatrix>);
    // compute S^-1. Assumes positive definite matrix
    void inverse();
    // compute S^-1/2. If an eigenvalue of S is smaller than thresh, the root will be discarded.
    void inverse_half(const double thresh = 1.0e-8);

    using Matrix_base<std::complex<double> >::copy_block;
    using Matrix_base<std::complex<double> >::get_block;
    void copy_block(const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const ZMatrix> o);
    // copy a block from a Matrix object to ZMatrix object and mutliply by coefficient coeff. If type = 0, coeff is purely real, else purely imaginary
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const double* data);
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::unique_ptr<double[]> o);
    void copy_real_block(const std::complex<double> a, const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const Matrix> o);

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

    // returns exp(*this)
    std::shared_ptr<ZMatrix> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<ZMatrix> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<ZMatrix> transpose() const;

    void zaxpy(const std::complex<double>, const ZMatrix&);
    void zaxpy(const std::complex<double>, const std::shared_ptr<const ZMatrix>);
    std::complex<double> zdotu(const ZMatrix&) const;
    std::complex<double> zdotu(const std::shared_ptr<const ZMatrix>) const;
    double norm() const;
    double rms() const;
    std::complex<double> trace() const;

    void zscal(const std::complex<double> a) { zscal_(size(), a, data(), 1); }
    void scale(const std::complex<double> a) { zscal(a); }

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

    std::complex<double> orthog(const std::list<std::shared_ptr<const ZMatrix> > o);

    // first parameter must be "R", "I", or "T" to print real part, imaginary part, or total
    void print(const std::string, const std::string in = "", const int size = 10) const;

    std::shared_ptr<ZMatrix> convert_real(const std::shared_ptr<const Matrix>);
};

}

#endif

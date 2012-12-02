//
// BAGEL - Parallel electron correlation program.
// Filename: matrix1e.h
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


#ifndef __src_util_zmatrix_h
#define __src_util_zmatrix_h

#include <cassert>
#include <string>
#include <algorithm>
#include <memory>
#include <list>
#include <src/util/f77.h>
#include <complex>

namespace bagel {

class ZMatrix { // Not to be confused with ZMatrix1e... at least for the moment
  protected:
    std::unique_ptr<std::complex<double>[]> data_;
    const int ndim_;
    const int mdim_;

  public:
    ZMatrix() : ndim_(0), mdim_(0) {};
    ZMatrix(const int n, const int m);
    ZMatrix(const ZMatrix&);
    ~ZMatrix();

    int size() const { return ndim_*mdim_; };
    int ndim() const { return ndim_; };
    int mdim() const { return mdim_; };
    std::complex<double>* data() const { return data_.get(); };
    std::complex<double>& data(const size_t i) { return *(data_.get()+i); };
    const std::complex<double>& data(const size_t i) const { return *(data_.get()+i); };
    std::complex<double>& element(int i, int j) { return *element_ptr(i, j); };
    std::complex<double>* element_ptr(int i, int j) { return data()+i+j*ndim_; };
    const std::complex<double>& element(int i, int j) const { return *element_ptr(i, j); };
    const std::complex<double>* element_ptr(int i, int j) const { return data()+i+j*ndim_; };

    void fill_upper();
    void symmetrize();
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
    void copy_block(const int nstart, const int mstart, const int ndim, const int mdim, const std::complex<double>* data);
    void copy_block(const int nstart, const int mstart, const int ndim, const int mdim, const std::unique_ptr<std::complex<double>[]> o);
    void copy_block(const int nstart, const int mstart, const int ndim, const int mdim, const std::shared_ptr<const ZMatrix> o);

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

    std::complex<double>& operator[](const size_t& i) { return data_[i]; };
    const std::complex<double>& operator[](const size_t& i) const { return data_[i]; };
    std::complex<double>& operator()(const size_t& i, const size_t& j) { return data_[i+j*ndim_]; };
    const std::complex<double>& operator()(const size_t& i, const size_t& j) const { return data_[i+j*ndim_]; };

    std::shared_ptr<ZMatrix> clone() const { return std::shared_ptr<ZMatrix>(new ZMatrix(ndim_, mdim_)); };

    // returns exp(*this)
    std::shared_ptr<ZMatrix> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<ZMatrix> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<ZMatrix> transpose() const;

    void zaxpy(const std::complex<double>, const ZMatrix&);
    void zaxpy(const std::complex<double>, const std::shared_ptr<const ZMatrix>);
    std::complex<double> zdotu(const ZMatrix&) const;
    std::complex<double> norm() const { return std::sqrt(zdotu(*this)); };
    std::complex<double> zdotu(const std::shared_ptr<const ZMatrix>) const;
    std::complex<double> rms() const;
    std::complex<double> trace() const;

    void zscal(const std::complex<double> a) { zscal_(size(), a, data(), 1); };
    void scale(const std::complex<double> a) { zscal(a); };

    void add_diag(const std::complex<double> a, const int i, const int j) { 
      assert(ndim_ == mdim_);
      for (int ii = i; ii != j; ++ii) element(ii,ii) += a;
    };
    void add_diag(const std::complex<double> a) { add_diag(a,0,ndim_); };
    // returns diagonal elements
    std::unique_ptr<std::complex<double>[]> diag() const;

    void fill(const std::complex<double> a) { std::fill(data(), data()+ndim_*mdim_, a); };
    void zero() { fill(0.0); };
    void unit() { fill(0.0); for (int i = 0; i != ndim_; ++i) element(i,i) = 1.0; assert(ndim_ == mdim_);};
    // purify a (near unitary) matrix to be unitary

    void purify_unitary();
    void purify_idempotent(const ZMatrix& s);
    void purify_redrotation(const int nclosed, const int nact, const int nvirt);

    std::complex<double> orthog(const std::list<std::shared_ptr<const ZMatrix> > o);

    void print(const std::string in = "", const int size = 10) const;
};

}

#endif

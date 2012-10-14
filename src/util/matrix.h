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


#ifndef __src_util_matrix_h
#define __src_util_matrix_h

#include <cassert>
#include <string>
#include <algorithm>
#include <memory>

namespace bagel {

class Matrix { // Not to be confused with Matrix... at least for the moment
  protected:
    std::unique_ptr<double[]> data_;
    int ndim_;
    int mdim_;

  public:
    Matrix() : ndim_(0), mdim_(0) {};
    Matrix(const int n, const int m);
    Matrix(const Matrix&);
    ~Matrix();

    int size() const { return ndim_*mdim_; };
    int ndim() const { return ndim_; };
    int mdim() const { return mdim_; };
    double* data() const { return data_.get(); };
    double& data(const size_t i) { return *(data_.get()+i); };
    const double& data(const size_t i) const { return *(data_.get()+i); };
    double& element(int i, int j) { return *element_ptr(i, j); };
    double* element_ptr(int i, int j) { return data()+i+j*ndim_; };
    const double& element(int i, int j) const { return *element_ptr(i, j); };
    const double* element_ptr(int i, int j) const { return data()+i+j*ndim_; };

    void fill_upper();
    void symmetrize();
    void antisymmetrize();
    std::shared_ptr<Matrix> resize(const int, const int) const;
    std::shared_ptr<Matrix> slice(const int, const int) const;
    std::shared_ptr<Matrix> expand() const;
    std::shared_ptr<Matrix> merge(const std::shared_ptr<const Matrix>) const;
    void diagonalize(double*);
    void svd(std::shared_ptr<Matrix>, std::shared_ptr<Matrix>);
    void inverse();

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

    std::shared_ptr<Matrix> clone() const { return std::shared_ptr<Matrix>(new Matrix(ndim_, mdim_)); };

    // returns exp(*this)
    std::shared_ptr<Matrix> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<Matrix> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<Matrix> transpose() const;

    void daxpy(const double, const Matrix&);
    void daxpy(const double, const std::shared_ptr<const Matrix>);
    double ddot(const Matrix&) const;
    double norm() const { return std::sqrt(ddot(*this)); };
    double ddot(const std::shared_ptr<const Matrix>) const;
    double rms() const;
    double trace() const;

    void dscal(const double a) { dscal_(size(), a, data(), 1); };
    void scale(const double a) { dscal(a); };

    void add_diag(const double a, const int i, const int j) { 
      assert(ndim_ == mdim_);
      for (int ii = i; ii != j; ++ii) element(ii,ii) += a;
    };
    void add_diag(const double a) { add_diag(a,0,ndim_); };
    // returns diagonal elements
    std::unique_ptr<double[]> diag() const;

    void fill(const double a) { std::fill(data(), data()+ndim_*mdim_, a); };
    void zero() { fill(0.0); };
    void unit() { fill(0.0); for (int i = 0; i != ndim_; ++i) element(i,i) = 1.0; assert(ndim_ == mdim_);};
    // purify a (near unitary) matrix to be unitary

    #if 0
    void purify_unitary();
    void purify_idempotent(const Matrix& s);
    void purify_redrotation(const int nclosed, const int nact, const int nvirt);

    double orthog(const std::list<std::shared_ptr<const Matrix> > o);
    #endif

    void print(const std::string in = "", const int size = 10) const;
};

}

#endif

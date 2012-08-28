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


#ifndef __src_scf_matrix1e_h
#define __src_scf_matrix1e_h

#include <cassert>
#include <src/scf/shell.h>
#include <src/scf/geometry.h>
#include <string>
#include <algorithm>
#include <memory>

namespace bagel {

class Matrix1e {
  protected:
    std::unique_ptr<double[]> data_;
    std::shared_ptr<const Geometry> geom_;
    int nbasis_;
    int ndim_;
    int mdim_;

    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int);
    virtual void init();

  public:
    Matrix1e() : nbasis_(0), ndim_(0), mdim_(0) {};
    Matrix1e(const std::shared_ptr<const Geometry>); 
    Matrix1e(const std::shared_ptr<const Geometry>, const int n, const int m);
    Matrix1e(const Matrix1e&); 
    ~Matrix1e();

    const std::shared_ptr<const Geometry> geom() const { return geom_; };

    int nbasis() const { return nbasis_; };
    int size() const { return nbasis_*nbasis_; };
    int ndim() const { return ndim_; }; 
    int mdim() const { return mdim_; }; 
    double* data() const { return data_.get(); };
    double& data(const size_t i) { return *(data_.get()+i); };
    const double& data(const size_t i) const { return *(data_.get()+i); };
    double& element(int i, int j) { return *element_ptr(i, j); };
    double* element_ptr(int i, int j) { return data()+i+j*nbasis_; };
    const double& element(int i, int j) const { return *element_ptr(i, j); };
    const double* element_ptr(int i, int j) const { return data()+i+j*nbasis_; };

    void fill_upper();
    void symmetrize();
    void antisymmetrize();
    std::shared_ptr<Matrix1e> resize(std::shared_ptr<const Geometry>, const int) const;
    std::shared_ptr<Matrix1e> slice(const int, const int) const;
    std::shared_ptr<Matrix1e> expand() const;
    std::shared_ptr<Matrix1e> merge(const std::shared_ptr<const Matrix1e>) const;
    void diagonalize(double*);
    void svd(std::shared_ptr<Matrix1e>, std::shared_ptr<Matrix1e>);
    void inverse();

    Matrix1e operator*(const Matrix1e&) const;
    Matrix1e& operator*=(const Matrix1e&);
    Matrix1e operator*(const double& a) const;
    Matrix1e operator/(const double& a) const;
    Matrix1e& operator*=(const double& a);
    Matrix1e& operator/=(const double& a);
    Matrix1e operator%(const Matrix1e&) const; // caution
    Matrix1e operator^(const Matrix1e&) const; // caution
    Matrix1e operator+(const Matrix1e&) const;
    Matrix1e& operator+=(const Matrix1e&);
    Matrix1e& operator-=(const Matrix1e&);
    Matrix1e& operator=(const Matrix1e&);
    Matrix1e operator-(const Matrix1e&) const;

    Matrix1e& operator/=(const Matrix1e&);
    Matrix1e operator/(const Matrix1e&) const;

    std::shared_ptr<Matrix1e> clone() const { return std::shared_ptr<Matrix1e>(new Matrix1e(geom_, ndim_, mdim_)); };

    // returns exp(*this)
    std::shared_ptr<Matrix1e> exp(const int deg = 6) const;
    // returns log(*this)
    std::shared_ptr<Matrix1e> log(const int deg = 6) const;
    // returns transpose(*this)
    std::shared_ptr<Matrix1e> transpose() const;

    void daxpy(const double, const Matrix1e&);
    void daxpy(const double, const std::shared_ptr<const Matrix1e>);
    double ddot(const Matrix1e&) const;
    double norm() const { return std::sqrt(ddot(*this)); };
    double ddot(const std::shared_ptr<const Matrix1e>) const;
    double rms() const;
    double trace() const;
    
    void dscal(const double a) { dscal_(size(), a, data(), 1); };
    void scale(const double a) { dscal(a); };

    void add_diag(const double a, const int i, const int j)
      { for (int ii = i; ii != j; ++ii) element(ii,ii) += a; };
    // returns diagonal elements
    std::unique_ptr<double[]> diag() const;

    void fill(const double a) { std::fill(data(), data()+nbasis_*nbasis_, a); };
    void zero() { fill(0.0); };
    void unit() { fill(0.0); for (int i = 0; i != ndim_; ++i) element(i,i) = 1.0; assert(ndim_ == mdim_);};
    // purify a (near unitary) matrix to be unitary
    void purify_unitary();
    void purify_idempotent(const Matrix1e& s);
    void purify_redrotation(const int nclosed, const int nact, const int nvirt);

    void print(const std::string in = "", const int size = 10) const;

    double orthog(const std::list<std::shared_ptr<const Matrix1e> > o);
};

}

#endif


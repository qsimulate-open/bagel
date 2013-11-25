//
// BAGEL - Parallel electron correlation program.
// Filename: matrix.cc
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


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/math/algo.h>
#include <src/math/matrix.h>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;


Matrix::Matrix(const int n, const int m , const bool loc) : Matrix_base<double>(n,m,loc) {
}


Matrix::Matrix(const Matrix& o) : Matrix_base<double>(o) {
}


Matrix::Matrix(Matrix&& o) : Matrix_base<double>(move(o)) {
}


#ifdef HAVE_SCALAPACK
Matrix::Matrix(const DistMatrix& o) : Matrix_base<double>(o.ndim(), o.mdim()) {
  setlocal_(o.local());
}
#endif


Matrix Matrix::operator+(const Matrix& o) const {
  Matrix out(*this);
  out.ax_plus_y(1.0, o);
  return out;
}


Matrix& Matrix::operator+=(const Matrix& o) {
  ax_plus_y(1.0, o);
  return *this;
}


Matrix& Matrix::operator-=(const Matrix& o) {
  ax_plus_y(-1.0, o);
  return *this;
}


Matrix& Matrix::operator=(const Matrix& o) {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  copy_n(o.data(), ndim_*mdim_, data());
  return *this;
}


Matrix& Matrix::operator=(Matrix&& o) {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  data_ = move(o.data_);
  return *this;
}


Matrix Matrix::operator-(const Matrix& o) const {
  Matrix out(*this);
  out.ax_plus_y(-1.0, o);
  return out;
}


Matrix Matrix::operator*(const Matrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.ndim());
  const int n = o.mdim();
  Matrix out(l, n, localized_);

#ifdef HAVE_SCALAPACK
  assert(localized_ == o.localized_);
  if (localized_) {
#endif
    dgemm_("N", "N", l, n, m, 1.0, data_, l, o.data_, o.ndim_, 0.0, out.data_, l);
#ifdef HAVE_SCALAPACK
  } else {
    unique_ptr<double[]> locala = getlocal();
    unique_ptr<double[]> localb = o.getlocal();
    unique_ptr<double[]> localc = out.getlocal();
    pdgemm_("N", "N", l, n, m, 1.0, locala.get(), desc_.get(), localb.get(), o.desc_.get(), 0.0, localc.get(), out.desc_.get());
    out.setlocal_(localc);
  }
#endif

  return out;
}


Matrix& Matrix::operator*=(const Matrix& o) {
  *this = *this * o;
  return *this;
}


Matrix Matrix::operator*(const double& a) const {
  Matrix out(*this);
  out *= a;
  return out;
}


Matrix Matrix::operator/(const double& a) const {
  Matrix out(*this);
  out /= a;
  return out;
}


Matrix& Matrix::operator*=(const double& a) {
  dscal_(ndim_*mdim_, a, data_, 1);
  return *this;
}
Matrix& Matrix::operator/=(const double& a) {
  *this *= 1.0/a;
  return *this;
}


Matrix Matrix::operator%(const Matrix& o) const {
  const int l = mdim_;
  const int m = ndim_;
  assert(ndim_ == o.ndim());
  const int n = o.mdim();
  Matrix out(l, n, localized_);

#ifdef HAVE_SCALAPACK
  assert(localized_ == o.localized_);
  if (localized_) {
#endif
    dgemm_("T", "N", l, n, m, 1.0, data_, m, o.data_, o.ndim_, 0.0, out.data_, l);
#ifdef HAVE_SCALAPACK
  } else {
    unique_ptr<double[]> locala = getlocal();
    unique_ptr<double[]> localb = o.getlocal();
    unique_ptr<double[]> localc = out.getlocal();
    pdgemm_("T", "N", l, n, m, 1.0, locala.get(), desc_.get(), localb.get(), o.desc_.get(), 0.0, localc.get(), out.desc_.get());
    out.setlocal_(localc);
  }
#endif

  return out;
}


Matrix Matrix::operator^(const Matrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.mdim());
  const int n = o.ndim();

  Matrix out(l, n, localized_);

#ifdef HAVE_SCALAPACK
  assert(localized_ == o.localized_);
  if (localized_) {
#endif
    dgemm_("N", "T", l, n, m, 1.0, data_, ndim_, o.data_, o.ndim_, 0.0, out.data_, l);
#ifdef HAVE_SCALAPACK
  } else {
    unique_ptr<double[]> locala = getlocal();
    unique_ptr<double[]> localb = o.getlocal();
    unique_ptr<double[]> localc = out.getlocal();
    pdgemm_("N", "T", l, n, m, 1.0, locala.get(), desc_.get(), localb.get(), o.desc_.get(), 0.0, localc.get(), out.desc_.get());
    out.setlocal_(localc);
  }
#endif

  return out;
}


Matrix Matrix::operator/(const Matrix& o) const {
  Matrix out(*this);
  out /= o;
  return out;
}


Matrix& Matrix::operator/=(const Matrix& o) {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);
  auto oiter = o.cbegin();
  for (auto& i : *this) {
    i /= *oiter++;
  }
  return *this;
}


void Matrix::diagonalize(double* eig) {
  if (ndim_ != mdim_) throw logic_error("illegal call of Matrix::diagonalize(double*)");
  const int n = ndim_;
  int info;

  // assume that the matrix is symmetric
  // the leading order (nbasis supplied)
#ifdef HAVE_SCALAPACK
  if (localized_ || n <= blocksize__) {
#endif
    unique_ptr<double[]> work(new double[n*6]);
    dsyev_("V", "L", n, data(), n, eig, work.get(), n*6, info);
    mpi__->broadcast(data(), n*n, 0);
#ifdef HAVE_SCALAPACK
  } else {
    const int localrow = get<0>(localsize_);
    const int localcol = get<1>(localsize_);

    unique_ptr<double[]> coeff(new double[localrow*localcol]);
    unique_ptr<double[]> local = getlocal();

    // first compute worksize
    double wsize;
    int liwork = 1;
    pdsyevd_("V", "U", n, local.get(), desc_.get(), eig, coeff.get(), desc_.get(), &wsize, -1, &liwork, 1, info);
    unique_ptr<int[]> iwork(new int[liwork]);
    wsize =  max(131072.0, wsize*2.0);

    const int lwork = round(wsize);
    unique_ptr<double[]> work(new double[lwork]);
    pdsyevd_("V", "U", n, local.get(), desc_.get(), eig, coeff.get(), desc_.get(), work.get(), lwork, iwork.get(), liwork, info);
    setlocal_(coeff);
  }
#endif

  if (info) throw runtime_error("dsyev/pdsyevd failed in Matrix");
}


void Matrix::svd(shared_ptr<Matrix> U, shared_ptr<Matrix> V) {
  assert(U->ndim() == ndim_ && U->mdim() == ndim_);
  assert(V->ndim() == mdim_ && V->mdim() == mdim_);
  const int lwork = 10*max(ndim_, mdim_);
  unique_ptr<double[]> work(new double[lwork]);
  unique_ptr<double[]> S(new double[min(ndim_, mdim_)]);
/*
  SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 $                   WORK, LWORK, INFO )
 */
  double* cblock = data();
  double* ublock = U->data();
  double* vblock = V->data();
  int info = 0;
  dgesvd_("A", "A", ndim_, mdim_, cblock, ndim_, S.get(), ublock, ndim_, vblock, mdim_, work.get(), lwork, info);
  if (info != 0) throw runtime_error("dgesvd failed in Matrix::svd");
}


shared_ptr<Matrix> Matrix::exp(const int deg) const {
  auto out = make_shared<Matrix>(ndim_, mdim_);
  Matrix buf(*this);
  assert(ndim_ == mdim_);

  for (int i = deg; i != 1; --i) {
    const double inv = 1.0/static_cast<double>(i);
    buf *= inv;
    for (int j = 0; j != ndim_; ++j) buf.element(j,j) += 1.0;
    *out = (*this)*buf;
    if (i != 1) buf = *out;
  }
  for (int j = 0; j != ndim_; ++j) out->element(j,j) += 1.0;
  return out;
}


shared_ptr<Matrix> Matrix::log(const int deg) const {
  auto out = make_shared<Matrix>(ndim_, mdim_);
  Matrix buf(*this);
  for (int j = 0; j != ndim_; ++j) buf.element(j,j) -= 1.0;
  assert(ndim_ == mdim_);

  for (int i = deg; i != 1; --i) {
    const double inv = -static_cast<double>(i-1)/static_cast<double>(i);
    buf *= inv;
    for (int j = 0; j != ndim_; ++j) buf.element(j,j) += 1.0;
    *out = (*this)*buf - buf;
    if (i != 1) buf = *out;
  }
  return out;
}


shared_ptr<Matrix> Matrix::transpose(const double factor) const {
  auto out = make_shared<Matrix>(mdim_, ndim_, localized_);
  mytranspose_(data_.get(), ndim_, mdim_, out->data(), factor);
  return out;
}


void Matrix::antisymmetrize() {
  assert(ndim_ == mdim_);
  shared_ptr<Matrix> trans = transpose();
  *this -= *trans;
  *this *= 0.5;
}


void Matrix::purify_unitary() {
  assert(ndim_ == mdim_);
  for (int i = 0; i != ndim_; ++i) {
    for (int j = 0; j != i; ++j) {
      const double a = ddot_(ndim_, &data_[i*ndim_], 1, &data_[j*ndim_], 1);
      daxpy_(ndim_, -a, &data_[j*ndim_], 1, &data_[i*ndim_], 1);
    }
    const double b = 1.0/std::sqrt(ddot_(ndim_, &data_[i*ndim_], 1, &data_[i*ndim_], 1));
    dscal_(ndim_, b, &data_[i*ndim_], 1);
  }
}


void Matrix::purify_redrotation(const int nclosed, const int nact, const int nvirt) {

#if 1
  for (int g = 0; g != nclosed; ++g)
    for (int h = 0; h != nclosed; ++h)
      element(h,g)=0.0;
  for (int g = 0; g != nact; ++g)
    for (int h = 0; h != nact; ++h)
      element(h+nclosed,g+nclosed)=0.0;
  for (int g = 0; g != nvirt; ++g)
    for (int h = 0; h != nvirt; ++h)
      element(h+nclosed+nact,g+nclosed+nact)=0.0;
  for (int i = 0; i != ndim_; ++i) {
    for (int j = 0; j != i; ++j) {
      const double ele = (element(j,i) - element(i,j)) * 0.5;
      element(j,i) = ele;
      element(i,j) = -ele;
    }
  }
#endif

}


void Matrix::purify_idempotent(const Matrix& s) {
  *this = *this * s * *this * 3.0 - *this * s * *this * s * *this * 2.0;
}


// in-place matrix inverse (practically we use buffer area)
void Matrix::inverse() {
  assert(ndim_ == mdim_);
  const int n = ndim_;
  shared_ptr<Matrix> buf = this->clone();
  buf->unit();

  int info;
  unique_ptr<int[]> ipiv(new int[n]);
  dgesv_(n, n, data(), n, ipiv.get(), buf->data(), n, info);
  if (info) throw runtime_error("dsysv failed in Matrix::inverse()");

  copy_n(buf->data(), n*n, data());
}


// this is a vector B, and solve AC = B, returns C
shared_ptr<Matrix> Matrix::solve(shared_ptr<const Matrix> A, const int n) const {
  Matrix As = *A;
  auto out = this->copy(); 
  assert(n <= out->ndim() && n <= A->ndim() && n <= A->mdim());

  unique_ptr<int[]> ipiv(new int[n]);
  int info;
  dgesv_(n, out->mdim(), As.data(), As.ndim(), ipiv.get(), out->data(), out->ndim(), info);
  if (info) throw runtime_error("DGESV failed");

  return out;
}


// compute S^{-1} using diagonalization
bool Matrix::inverse_symmetric(const double thresh) {
  assert(ndim_ == mdim_);
  const int n = ndim_;
  unique_ptr<double[]> vec(new double[n]);
  diagonalize(vec.get());

  for (int i = 0; i != n; ++i) {
    double s = vec[i] > thresh ? 1.0/std::sqrt(vec[i]) : 0.0;
    dscal_(n, s, data_.get()+i*n, 1);
  }
  *this = *this ^ *this;
  vector<double> rm;
  for (int i = 0; i != n; ++i)
    if (vec[i] < thresh) rm.push_back(vec[i]);
#ifndef NDEBUG
  if (!rm.empty())
    cout << "    - linear dependency detected: " << setw(4) << rm.size() << " / " << setw(4) << n <<
            "    min eigenvalue: " << setw(14) << scientific << setprecision(4) << *min_element(rm.begin(), rm.end()) <<
            "    max eigenvalue: " << setw(14) << scientific << setprecision(4) << *max_element(rm.begin(), rm.end()) << fixed << endl;
#endif
  return rm.empty();
}


// compute S^{-1/2}
bool Matrix::inverse_half(const double thresh) {
  assert(ndim_ == mdim_);
  const int n = ndim_;
  unique_ptr<double[]> vec(new double[n]);

#ifdef HAVE_SCALAPACK
  if (localized_) {
#endif
    diagonalize(vec.get());
    for (int i = 0; i != n; ++i) {
      double s = vec[i] > thresh ? 1.0/std::sqrt(std::sqrt(vec[i])) : 0.0;
      dscal_(n, s, data_.get()+i*n, 1);
    }
    *this = *this ^ *this;
#ifdef HAVE_SCALAPACK
  } else {
    unique_ptr<double[]> scal(new double[n]);
    shared_ptr<DistMatrix> dist = distmatrix();
    dist->diagonalize(vec.get());
    for (int i = 0; i != n; ++i)
      scal[i] = vec[i] > thresh ? 1.0/std::sqrt(std::sqrt(vec[i])) : 0.0;
    dist->scale(scal.get());
    *this = *(*dist ^ *dist).matrix();
  }
#endif

  vector<double> rm;
  for (int i = 0; i != n; ++i)
    if (vec[i] < thresh) rm.push_back(vec[i]);
  if (!rm.empty())
    cout << "    - linear dependency detected: " << setw(4) << rm.size() << " / " << setw(4) << n <<
            "    min eigenvalue: " << setw(14) << scientific << setprecision(4) << *min_element(rm.begin(), rm.end()) <<
            "    max eigenvalue: " << setw(14) << scientific << setprecision(4) << *max_element(rm.begin(), rm.end()) << fixed << endl;
  return rm.empty();
}

// compute Hermitian square root, S^{1/2}
void Matrix::sqrt() {
  assert(ndim_ == mdim_);
  const int n = ndim_;
  unique_ptr<double[]> vec(new double[n]);
  diagonalize(vec.get());

  for (int i = 0; i != n; ++i) {
    if (vec[i] < 0.0) throw runtime_error("Matrix::sqrt() called, but this matrix is not positive definite");
    double s = std::sqrt(std::sqrt(vec[i]));
    dscal_(n, s, data_.get()+i*n, 1);
  }

  *this = *this ^ *this;
}


void Matrix::print(const string name, const size_t size) const {

  cout << "++++ " + name + " ++++" << endl;
  for (int i = 0; i != min(size,ndim_); ++i) {
    for (int j = 0; j != min(size,mdim_); ++j) {
      cout << fixed << setw(12) << setprecision(9) << data_[j * ndim_ + i]  << " ";
    }
    cout << endl;
  }

}


#ifdef HAVE_SCALAPACK
shared_ptr<DistMatrix> Matrix::distmatrix() const {
  return make_shared<DistMatrix>(*this);
}
#else
shared_ptr<const Matrix> Matrix::distmatrix() const {
  shared_ptr<const Matrix> out = shared_from_this();
  return out;
}
#endif


#ifndef HAVE_SCALAPACK
shared_ptr<const Matrix> Matrix::form_density_rhf(const int n, const int offset) const {
  shared_ptr<const Matrix> tmp = this->slice(offset, offset+n);
  auto out = make_shared<Matrix>(*tmp ^ *tmp);
  *out *= 2.0;
  return out;
}
#endif

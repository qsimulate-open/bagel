//
// BAGEL - Parallel electron correlation program.
// Filename: matrix1e.cc
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


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/util/matrix.h>
#include <cassert>
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace bagel;

Matrix::Matrix(const int n, const int m)
 : data_(new double[n*m]) {
  ndim_ = n;
  mdim_ = m;
  zero();
}


Matrix::Matrix(const Matrix& o)
 : data_(new double[o.ndim_*o.mdim_]), ndim_(o.ndim_), mdim_(o.mdim_) {
  copy(o.data(), o.data() + ndim_*mdim_, data());
}

Matrix::~Matrix() {
}


shared_ptr<Matrix> Matrix::resize(const int n, const int m) const {
  shared_ptr<Matrix> out(new Matrix(n, m));
  for (int i = 0; i != mdim_; ++i) {
    for (int j = 0; j != ndim_; ++j) {
      out->data_[j+i*n] = data_[j+i*ndim_];
    }
  }
  return out;
}


shared_ptr<Matrix> Matrix::slice(const int start, const int fence) const {
  shared_ptr<Matrix> out(new Matrix(ndim_, fence - start));
  assert(fence <= ndim_);

  copy(data_.get()+start*ndim_, data_.get()+fence*ndim_, out->data_.get());
  return out;
}


shared_ptr<Matrix> Matrix::merge(const shared_ptr<const Matrix> o) const {
  assert(ndim_ == o->ndim_);

  shared_ptr<Matrix> out(new Matrix(ndim_, mdim_ + o->mdim_));

  copy(data_.get(), data_.get() + ndim_*mdim_, out->data_.get());
  copy(o->data_.get(), o->data_.get()+o->ndim_*o->mdim_, out->data_.get()+ndim_*mdim_);
  return out;
}


Matrix Matrix::operator+(const Matrix& o) const {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);

  Matrix out(*this);
  const int size = ndim_ * mdim_;
  const double* odata = o.data();
  double* outdata = out.data();

  daxpy_(size, 1.0, data(), 1, outdata, 1);
  daxpy_(size, 1.0, odata, 1, outdata, 1);

  return out;
}


Matrix& Matrix::operator+=(const Matrix& o) {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);

  daxpy_(ndim_*mdim_, 1.0, o.data(), 1, data(), 1);
  return *this;
}


Matrix& Matrix::operator-=(const Matrix& o) {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);

  daxpy_(ndim_*mdim_, -1.0, o.data(), 1, data(), 1);
  return *this;
}


Matrix& Matrix::operator=(const Matrix& o) {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  dcopy_(ndim_*mdim_, o.data(), 1, data(), 1);
  return *this;
}


Matrix Matrix::operator-(const Matrix& o) const {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);
  Matrix out(*this);

  const int size = ndim_ * mdim_;

  const double* odata = o.data();
  double* outdata = out.data();
  daxpy_(size, 1.0, data(), 1, outdata, 1);
  daxpy_(size, -1.0, odata, 1, outdata, 1);

  return out;
}


Matrix Matrix::operator*(const Matrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.ndim());
  const int n = o.mdim();

  Matrix out(l, n);
  const double* odata = o.data();
  double* outdata = out.data();

  dgemm_("N", "N", l, n, m, 1.0, data(), l, odata, o.ndim_, 0.0, outdata, l);

  return out;
}


Matrix& Matrix::operator*=(const Matrix& o) {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.ndim());
  const int n = o.mdim();

  Matrix out(l, n);
  const double* odata = o.data();
  double* outdata = out.data();

  dgemm_("N", "N", l, n, m, 1.0, data(), l, odata, o.ndim_, 0.0, outdata, l);

  *this = out;
  return *this;
}


Matrix Matrix::operator*(const double& a) const {
  Matrix out(*this);
  dscal_(ndim_*mdim_, a, out.data(), 1);
  return out;
}


Matrix Matrix::operator/(const double& a) const {
  Matrix out(*this);
  dscal_(ndim_*mdim_, 1.0/a, out.data(), 1);
  return out;
}


Matrix& Matrix::operator*=(const double& a) {
  dscal_(ndim_*mdim_, a, data_, 1);
  return *this;
}
Matrix& Matrix::operator/=(const double& a) {
  dscal_(ndim_*mdim_, 1.0/a, data_, 1);
  return *this;
}


Matrix Matrix::operator%(const Matrix& o) const {
  const int l = mdim_;
  const int m = ndim_;
  assert(ndim_ == o.ndim());
  const int n = o.mdim();

  Matrix out(l, n);

  const double* odata = o.data();
  double* outdata = out.data();

  dgemm_("T", "N", l, n, m, 1.0, data(), m, odata, o.ndim_, 0.0, outdata, l);

  return out;
}


Matrix Matrix::operator^(const Matrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  const int n = o.ndim();

  Matrix out(l, n);

  assert(mdim_ == o.mdim());
  const double* odata = o.data();
  double* outdata = out.data();

  dgemm_("N", "T", l, n, m, 1.0, data(), ndim_, odata, o.ndim_, 0.0, outdata, l);

  return out;
}


Matrix Matrix::operator/(const Matrix& o) const {
  Matrix out(*this);
  out /= o;
  return out;
}


Matrix& Matrix::operator/=(const Matrix& o) {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);
  for (int i = 0; i != size(); ++i) {
    data(i) /= o.data(i);
  }
  return *this;
}


void Matrix::diagonalize(double* eig) {
  if (ndim_ != mdim_) throw logic_error("illegal call of Matrix::diagonalize(double*)"); 
  const int n = ndim_;
  // assume that the matrix is symmetric
  // the leading order (nbasis supplied)
  int info;
  unique_ptr<double[]> work(new double[n*6]);
  dsyev_("V", "L", n, data(), n, eig, work.get(), n*6, info);
  if(info) throw runtime_error("diagonalize failed");
}


#if 0 // All the following may not be needed just yet
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
  dgesvd_("A", "A", &ndim_, &mdim_, cblock, &ndim_, S.get(), ublock, &ndim_, vblock, &mdim_, work.get(), &lwork, &info);
  if (info != 0) throw runtime_error("dgesvd failed in Matrix::svd");
}
#endif


void Matrix::daxpy(const double a, const Matrix& o) {
  daxpy_(ndim_*mdim_, a, o.data(), 1, data(), 1);
}


void Matrix::daxpy(const double a, const std::shared_ptr<const Matrix> o) {
  daxpy_(ndim_*mdim_, a, o->data(), 1, data(), 1);
}


double Matrix::ddot(const Matrix& o) const {
  return ddot_(ndim_*mdim_, data(), 1, o.data(), 1);
}


double Matrix::ddot(const std::shared_ptr<const Matrix> o) const {
  return ddot_(ndim_*mdim_, data(), 1, o->data(), 1);
}


double Matrix::rms() const {
  return ::sqrt(ddot(*this) / (ndim_ * mdim_));
}


double Matrix::trace() const {
  double out = 0.0;
  assert(ndim_ == mdim_);
  for (int i = 0; i != ndim_; ++i) out += data_[i * ndim_ + i];
  return out;
}


#if 0 // Only makes sense for square matrices
shared_ptr<Matrix> Matrix::exp(const int deg) const {
  shared_ptr<Matrix> out(new Matrix(geom_, ndim_, mdim_));
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
  shared_ptr<Matrix> out(new Matrix(geom_, ndim_, mdim_));
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
#endif


unique_ptr<double[]> Matrix::diag() const {
  if (ndim_ != mdim_) throw logic_error("illegal call of Matrix::diag()");
  unique_ptr<double[]> out(new double[ndim_]);
  for (int i = 0; i != ndim_; ++i) {
    out[i] = element(i,i);
  }
  return move(out);
}


shared_ptr<Matrix> Matrix::transpose() const {
  shared_ptr<Matrix> out(new Matrix(mdim_, ndim_));

#if 0
  double *data = data_.get();
  double *odata = out->data();
  for(int i = 0; i < ndim_; ++i) {
    for(int j = 0; j < mdim_; ++j) {
      odata[i*mdim_ + j] = data[i + j*ndim_];
    }
  }
#else
  mytranspose_(data_.get(), &ndim_, &mdim_, out->data()); 
#endif

  return out;
}


void Matrix::symmetrize() {
  assert(ndim_ == mdim_);
  const int n = mdim_;
  for (int i = 0; i != n; ++i)
    for (int j = i+1; j != n; ++j)
      data_[i+j*n] = data_[j+i*n] = 0.5*(data_[i+j*n]+data_[j+i*n]);
}


void Matrix::antisymmetrize() {
  assert(ndim_ == mdim_);

  shared_ptr<Matrix> trans = transpose();
  *this -= *trans;
  *this *= 0.5;
}


#if 0
void Matrix::purify_unitary() {
  assert(ndim_ == mdim_);
#if 0
  Matrix buf(*this ^ *this);
  const int lwork = 5*ndim_;
  unique_ptr<double[]> work(new double[lwork]);
  unique_ptr<double[]> vec(new double[ndim_]);
  int info;
  dsyev_("V", "U", ndim_, buf.data(), ndim_, vec.get(), work.get(), lwork, info);
  if (info) throw runtime_error("dsyev failed in Matrix::purify_unitary");
  if (vec[0] < 0.95)        cout << "   --- smallest eigenvalue in purify_unitary() " << vec[0] << endl;
  if (vec[ndim_-1] > 1.05)  cout << "   --- largest eigenvalue in purify_unitary() " << vec[ndim_-1] << endl;
  for (int i = 0; i != ndim_; ++i) {
    for (int j = 0; j != ndim_; ++j) {
      buf.element(j,i) /= std::sqrt(std::sqrt(vec[i]));
    }
  }
  *this = ((buf ^ buf) * *this);
#else
  // Schmidt orthogonalization
  for (int i = 0; i != ndim_; ++i) {
    for (int j = 0; j != i; ++j) {
      const double a = ddot_(ndim_, &data_[i*ndim_], 1, &data_[j*ndim_], 1);
      daxpy_(ndim_, -a, &data_[j*ndim_], 1, &data_[i*ndim_], 1);
    }
    const double b = 1.0/sqrt(ddot_(ndim_, &data_[i*ndim_], 1, &data_[i*ndim_], 1));
    dscal_(ndim_, b, &data_[i*ndim_], 1);
  }
#endif

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
  for (int i = 0; i != nbasis_; ++i) {
    for (int j = 0; j != i; ++j) {
      const double ele = (element(j,i) - element(i,j)) * 0.5;
      element(j,i) = ele;
      element(i,j) = -ele;
    }
  }
#endif

}
#endif


#if 0 // not sure I need all of this
void Matrix::purify_idempotent(const Matrix& s) {
  *this = *this * s * *this * 3.0 - *this * s * *this * s * *this * 2.0;
}


double Matrix::orthog(const std::list<std::shared_ptr<const Matrix> > o) {
  for (auto& it : o) {
    const double m = this->ddot(it);
    this->daxpy(-m, it);
  }
  const double n = norm();
  *this /= n;
  return n;
}
#endif


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


void Matrix::print(const string name, const int size) const {

  cout << "++++ " + name + " ++++" << endl;
  for (int i = 0; i != min(size,ndim_); ++i) {
    for (int j = 0; j != min(size,mdim_); ++j) {
      cout << fixed << setw(12) << setprecision(9) << data_[j * ndim_ + i]  << " ";
    }
    cout << endl;
  }

}

void Matrix::copy_block(const int ndim_i, const int mdim_i, const int ndim, const int mdim, const double* data) {
  for (int i = mdim_i, j = 0; i != mdim_i + mdim ; ++i, ++j)
    copy_n(data + j*ndim, ndim, data_.get() + ndim_i + i*ndim_);
}


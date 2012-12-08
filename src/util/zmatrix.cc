//
// BAGEL - Parallel electron correlation program.
// Filename: zmatrix.cc
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


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/util/zmatrix.h>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

ZMatrix::ZMatrix(const int n, const int m) : data_(new complex<double>[n*m]), ndim_(n), mdim_(m) {
  zero();
}


ZMatrix::ZMatrix(const ZMatrix& o) : data_(new complex<double>[o.ndim_*o.mdim_]), ndim_(o.ndim_), mdim_(o.mdim_) {
  copy_n(o.data(), ndim_*mdim_, data());
}

ZMatrix::~ZMatrix() {
}


shared_ptr<ZMatrix> ZMatrix::cut(const int n) const {
  assert(n <= ndim_);
  shared_ptr<ZMatrix> out(new ZMatrix(n, mdim_));
  for (int i = 0; i != mdim_; ++i)
    for (int j = 0; j != n; ++j)
      out->data_[j+i*n] = data_[j+i*ndim_];
  return out;
}



shared_ptr<ZMatrix> ZMatrix::resize(const int n, const int m) const {
  assert(n >= ndim_ && m >= mdim_);
  shared_ptr<ZMatrix> out(new ZMatrix(n, m));
  for (int i = 0; i != mdim_; ++i) {
    for (int j = 0; j != ndim_; ++j) {
      out->data_[j+i*n] = data_[j+i*ndim_];
    }
  }
  return out;
}


shared_ptr<ZMatrix> ZMatrix::slice(const int start, const int fence) const {
  shared_ptr<ZMatrix> out(new ZMatrix(ndim_, fence - start));
  assert(fence <= ndim_);

  copy(data_.get()+start*ndim_, data_.get()+fence*ndim_, out->data_.get());
  return out;
}


shared_ptr<ZMatrix> ZMatrix::merge(const shared_ptr<const ZMatrix> o) const {
  assert(ndim_ == o->ndim_);

  shared_ptr<ZMatrix> out(new ZMatrix(ndim_, mdim_ + o->mdim_));

  copy(data_.get(), data_.get() + ndim_*mdim_, out->data_.get());
  copy(o->data_.get(), o->data_.get()+o->ndim_*o->mdim_, out->data_.get()+ndim_*mdim_);
  return out;
}


ZMatrix ZMatrix::operator+(const ZMatrix& o) const {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);

  ZMatrix out(*this);
  const int size = ndim_ * mdim_;
  const complex<double>* odata = o.data();
  complex<double>* outdata = out.data();

  zaxpy_(size, 1.0, odata, 1, outdata, 1);

  return out;
}


ZMatrix& ZMatrix::operator+=(const ZMatrix& o) {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);

  zaxpy_(ndim_*mdim_, 1.0, o.data(), 1, data(), 1);
  return *this;
}


ZMatrix& ZMatrix::operator-=(const ZMatrix& o) {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);

  zaxpy_(ndim_*mdim_, -1.0, o.data(), 1, data(), 1);
  return *this;
}


ZMatrix& ZMatrix::operator=(const ZMatrix& o) {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  copy_n(o.data(), ndim_*mdim_, data());
  return *this;
}


ZMatrix ZMatrix::operator-(const ZMatrix& o) const {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);
  ZMatrix out(*this);

  const int size = ndim_ * mdim_;

  const complex<double>* odata = o.data();
  complex<double>* outdata = out.data();
  zaxpy_(size, -1.0, odata, 1, outdata, 1);

  return out;
}


ZMatrix ZMatrix::operator*(const ZMatrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.ndim());
  const int n = o.mdim();

  ZMatrix out(l, n);
  const complex<double>* odata = o.data();
  complex<double>* outdata = out.data();

  zgemm3m_("N", "N", l, n, m, 1.0, data(), l, odata, o.ndim_, 0.0, outdata, l);

  return out;
}


ZMatrix& ZMatrix::operator*=(const ZMatrix& o) {
  *this = *this * o;
  return *this;
}


ZMatrix ZMatrix::operator*(const complex<double>& a) const {
  ZMatrix out(*this);
  zscal_(ndim_*mdim_, a, out.data(), 1);
  return out;
}


ZMatrix ZMatrix::operator/(const complex<double>& a) const {
  ZMatrix out(*this);
  zscal_(ndim_*mdim_, 1.0/a, out.data(), 1);
  return out;
}


ZMatrix& ZMatrix::operator*=(const complex<double>& a) {
  zscal_(ndim_*mdim_, a, data_, 1);
  return *this;
}
ZMatrix& ZMatrix::operator/=(const complex<double>& a) {
  zscal_(ndim_*mdim_, 1.0/a, data_, 1);
  return *this;
}


ZMatrix ZMatrix::operator%(const ZMatrix& o) const {
  const int l = mdim_;
  const int m = ndim_;
  assert(ndim_ == o.ndim());
  const int n = o.mdim();

  ZMatrix out(l, n);

  const complex<double>* odata = o.data();
  complex<double>* outdata = out.data();

  zgemm3m_("H", "N", l, n, m, 1.0, data(), m, odata, o.ndim_, 0.0, outdata, l);

  return out;
}


ZMatrix ZMatrix::operator^(const ZMatrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  const int n = o.ndim();

  ZMatrix out(l, n);

  assert(mdim_ == o.mdim());
  const complex<double>* odata = o.data();
  complex<double>* outdata = out.data();

  zgemm3m_("N", "H", l, n, m, 1.0, data(), ndim_, odata, o.ndim_, 0.0, outdata, l);

  return out;
}


ZMatrix ZMatrix::operator/(const ZMatrix& o) const {
  ZMatrix out(*this);
  out /= o;
  return out;
}


ZMatrix& ZMatrix::operator/=(const ZMatrix& o) {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);
  for (int i = 0; i != size(); ++i) {
    data(i) /= o.data(i);
  }
  return *this;
}


void ZMatrix::diagonalize(double* eig) {
  if (ndim_ != mdim_) throw logic_error("illegal call of ZMatrix::diagonalize(complex<double>*)"); 
  const int n = ndim_;
  // assume that the matrix is symmetric
  // the leading order (nbasis supplied)
  int info;
  unique_ptr<double[]> rwork(new double[max(1, 3*ndim_-2)]);

  unique_ptr<complex<double>[]> work(new complex<double>[n*6]);
  zheev_("V", "L", n, data(), n, eig, work.get(), n*6, rwork.get(), info);
  if(info) throw runtime_error("diagonalize failed");

  // in order to be consistent with signs
  mpi__->broadcast(data(), n*n, 0);
}


void ZMatrix::svd(shared_ptr<ZMatrix> U, shared_ptr<ZMatrix> V) {
  assert(U->ndim() == ndim_ && U->mdim() == ndim_);
  assert(V->ndim() == mdim_ && V->mdim() == mdim_);
  const int lwork = 10*max(ndim_, mdim_);
  unique_ptr<double[]> rwork(new double[5*max(ndim_, mdim_)]);
  unique_ptr<complex<double>[]> work(new complex<double>[lwork]);
  unique_ptr<double[]> S(new double[min(ndim_, mdim_)]);
/*
  SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 $                   WORK, LWORK, RWORK, INFO )
 */
  complex<double>* cblock = data();
  complex<double>* ublock = U->data();
  complex<double>* vblock = V->data();
  int info = 0;
  zgesvd_("A", "A", &ndim_, &mdim_, cblock, &ndim_, S.get(), ublock, &ndim_, vblock, &mdim_, work.get(), &lwork, rwork.get(), &info);
  if (info != 0) throw runtime_error("zgesvd failed in ZMatrix::svd");
}


void ZMatrix::zaxpy(const complex<double> a, const ZMatrix& o) {
  zaxpy_(ndim_*mdim_, a, o.data(), 1, data(), 1);
}


void ZMatrix::zaxpy(const complex<double> a, const std::shared_ptr<const ZMatrix> o) {
  zaxpy_(ndim_*mdim_, a, o->data(), 1, data(), 1);
}


complex<double> ZMatrix::zdotu(const ZMatrix& o) const {
  return zdotu_(ndim_*mdim_, data(), 1, o.data(), 1);
}


complex<double> ZMatrix::zdotu(const std::shared_ptr<const ZMatrix> o) const {
  return zdotu_(ndim_*mdim_, data(), 1, o->data(), 1);
}


double ZMatrix::norm() const {
  complex<double> n = zdotu(*this);
  assert(fabs(n.imag()) < 1.0e-10);
  return n.real();
}


double ZMatrix::rms() const {
  return norm()/std::sqrt(ndim_ * mdim_);
}


complex<double> ZMatrix::trace() const {
  complex<double> out = 0.0;
  assert(ndim_ == mdim_);
  for (int i = 0; i != ndim_; ++i) out += data_[i * ndim_ + i];
  return out;
}


shared_ptr<ZMatrix> ZMatrix::exp(const int deg) const {
  shared_ptr<ZMatrix> out(new ZMatrix(ndim_, mdim_));
  ZMatrix buf(*this);
  assert(ndim_ == mdim_);

  for (int i = deg; i != 1; --i) {
    const complex<double> inv = 1.0/static_cast<complex<double>>(i);
    buf *= inv;
    for (int j = 0; j != ndim_; ++j) buf.element(j,j) += 1.0;
    *out = (*this)*buf;
    if (i != 1) buf = *out;
  }
  for (int j = 0; j != ndim_; ++j) out->element(j,j) += 1.0;
  return out;
}


shared_ptr<ZMatrix> ZMatrix::log(const int deg) const {
  shared_ptr<ZMatrix> out(new ZMatrix(ndim_, mdim_));
  ZMatrix buf(*this);
  for (int j = 0; j != ndim_; ++j) buf.element(j,j) -= 1.0;
  assert(ndim_ == mdim_);

  for (int i = deg; i != 1; --i) {
    const complex<double> inv = -static_cast<complex<double>>(i-1)/static_cast<complex<double>>(i);
    buf *= inv;
    for (int j = 0; j != ndim_; ++j) buf.element(j,j) += 1.0;
    *out = (*this)*buf - buf;
    if (i != 1) buf = *out;
  }
  return out;
}


unique_ptr<complex<double>[]> ZMatrix::diag() const {
  if (ndim_ != mdim_) throw logic_error("illegal call of ZMatrix::diag()");
  unique_ptr<complex<double>[]> out(new complex<double>[ndim_]);
  for (int i = 0; i != ndim_; ++i) {
    out[i] = element(i,i);
  }
  return move(out);
}


shared_ptr<ZMatrix> ZMatrix::transpose() const {
  shared_ptr<ZMatrix> out(new ZMatrix(mdim_, ndim_));

#if 1
  complex<double> *data = data_.get();
  complex<double> *odata = out->data();
  for(int i = 0; i < ndim_; ++i) {
    for(int j = 0; j < mdim_; ++j) {
      odata[i*mdim_ + j] = data[i + j*ndim_];
    }
  }
#else
//  mytranspose_(data_.get(), &ndim_, &mdim_, out->data()); 
#endif

  return out;
}


void ZMatrix::fill_upper() {
  assert(ndim_ == mdim_);
  for (int i = 0; i != mdim_; ++i)
    for (int j = i+1; j != ndim_; ++j)
      data_[i+j*ndim_] = data_[j+i*ndim_];
}


void ZMatrix::symmetrize() {
  assert(ndim_ == mdim_);
  const int n = mdim_;
  for (int i = 0; i != n; ++i)
    for (int j = i+1; j != n; ++j)
      data_[i+j*n] = data_[j+i*n] = 0.5*(data_[i+j*n]+data_[j+i*n]);
}


void ZMatrix::antisymmetrize() {
  assert(ndim_ == mdim_);

  shared_ptr<ZMatrix> trans = transpose();
  *this -= *trans;
  *this *= 0.5;
}


void ZMatrix::purify_unitary() {
  assert(ndim_ == mdim_);
#if 0
  ZMatrix buf(*this ^ *this);
  const int lwork = 5*ndim_;
  unique_ptr<complex<double>[]> work(new complex<double>[lwork]);
  unique_ptr<complex<double>[]> vec(new complex<double>[ndim_]);
  int info;
  dsyev_("V", "U", ndim_, buf.data(), ndim_, vec.get(), work.get(), lwork, info);
  if (info) throw runtime_error("dsyev failed in ZMatrix::purify_unitary");
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
      const complex<double> a = zdotu_(ndim_, &data_[i*ndim_], 1, &data_[j*ndim_], 1);
      zaxpy_(ndim_, -a, &data_[j*ndim_], 1, &data_[i*ndim_], 1);
    }
    const complex<double> b = 1.0/sqrt(zdotu_(ndim_, &data_[i*ndim_], 1, &data_[i*ndim_], 1));
    zscal_(ndim_, b, &data_[i*ndim_], 1);
  }
#endif

}


void ZMatrix::purify_redrotation(const int nclosed, const int nact, const int nvirt) {

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
      const complex<double> ele = (element(j,i) - element(i,j)) * 0.5;
      element(j,i) = ele;
      element(i,j) = -ele;
    }
  }
#endif

}


void ZMatrix::purify_idempotent(const ZMatrix& s) {
  *this = *this * s * *this * 3.0 - *this * s * *this * s * *this * 2.0;
}


complex<double> ZMatrix::orthog(const std::list<std::shared_ptr<const ZMatrix> > o) {
  for (auto& it : o) {
    const complex<double> m = this->zdotu(it);
    this->zaxpy(-m, it);
  }
  const complex<double> n = norm();
  *this /= n;
  return n;
}


// in-place matrix inverse (practically we use buffer area)
void ZMatrix::inverse() {
  assert(ndim_ == mdim_);
  const int n = ndim_;
  shared_ptr<ZMatrix> buf = this->clone();
  buf->unit();

  int info;
  unique_ptr<int[]> ipiv(new int[n]);
  zgesv_(n, n, data(), n, ipiv.get(), buf->data(), n, info);
  if (info) throw runtime_error("dsysv failed in ZMatrix::inverse()");

  copy_n(buf->data(), n*n, data());
}


// compute S^{-1/2}
void ZMatrix::inverse_half(const double thresh) {
  assert(ndim_ == mdim_);
  const int n = ndim_;
  unique_ptr<double[]> vec(new double[n]);
  diagonalize(vec.get());

  for (int i = 0; i != n; ++i) {
    double s = vec[i] > thresh ? 1.0/sqrt(sqrt(vec[i])) : 0.0;
    zscal_(n, s, data_.get()+i*n, 1);
  }

#ifndef NDEBUG
  for (int i = 0; i != n; ++i)
    if (vec[i] < thresh) cout << " throwing out " << setprecision(20) << vec[i] << endl;
#endif

  *this = *this ^ *this;

}


void ZMatrix::print(const string name, const int size) const {

  cout << "++++ " + name + " ++++" << endl;
  for (int i = 0; i != min(size,ndim_); ++i) {
    for (int j = 0; j != min(size,mdim_); ++j) {
      cout << fixed << setw(12) << setprecision(9) << data_[j * ndim_ + i]  << " ";
    }
    cout << endl;
  }

}


void ZMatrix::copy_block(const int ndim_i, const int mdim_i, const int ndim, const int mdim, const complex<double>* data) {
  for (int i = mdim_i, j = 0; i != mdim_i + mdim ; ++i, ++j)
    copy_n(data + j*ndim, ndim, data_.get() + ndim_i + i*ndim_);
}
 

void ZMatrix::copy_block(const int ndim_i, const int mdim_i, const int ndim, const int mdim, const unique_ptr<complex<double>[]> data) { copy_block(ndim_i, mdim_i, ndim, mdim, data.get()); }
void ZMatrix::copy_block(const int ndim_i, const int mdim_i, const int ndim, const int mdim, const shared_ptr<const ZMatrix> data) { copy_block(ndim_i, mdim_i, ndim, mdim, data->data()); }

void ZMatrix::copy_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const double* data) {
//if (type == 0) {
    for (int i = mdim_i, j = 0; i != mdim_i + mdim ; ++i, ++j) {
      for (int k = ndim_i, l = 0; k != ndim_i + ndim ; ++k, ++l) {
        element(k,i) = a * *(data + j*ndim + l);
//      element(k,i) = a * (*(data + j*ndim + l), 0);
      } 
    } 
#if 0
  } else {
    for (int i = mdim_i, j = 0; i != mdim_i + mdim ; ++i, ++j) {
      for (int k = ndim_i, l = 0; k != ndim_i + ndim ; ++k, ++l) {
        element(k,i) = (0, coeff * *(data + j*ndim + l));
      } 
    } 
  }
#endif
} 
 
void ZMatrix::copy_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const unique_ptr<double[]> data) { copy_real_block(a, ndim_i, mdim_i, ndim, mdim, data.get()); }

void ZMatrix::copy_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const shared_ptr<const Matrix> data) {
  assert(ndim == data->ndim() && mdim == data->mdim());
  copy_real_block(a, ndim_i, mdim_i, ndim, mdim, data->data());
}

shared_ptr<ZMatrix> ZMatrix::convert_real(const shared_ptr<const Matrix> in) {
  shared_ptr<ZMatrix> out(new ZMatrix(in->ndim(), in->mdim()));
  for (int i = 0; i != in->size(); ++i) {
    out->data_[i] = (in->data(i), 0);
  }
  return out;
}

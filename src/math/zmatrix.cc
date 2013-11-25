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


#include <iostream>
#include <iomanip>
#include <cmath>
#include <src/math/zmatrix.h>

using namespace std;
using namespace bagel;

ZMatrix::ZMatrix(const int n, const int m, const bool loc) : Matrix_base<complex<double>>(n, m, loc) {
}


ZMatrix::ZMatrix(const ZMatrix& o) : Matrix_base<complex<double>>(o) {
}


ZMatrix::ZMatrix(ZMatrix&& o) : Matrix_base<complex<double>>(move(o)) {
}


ZMatrix::ZMatrix(const Matrix& r, const Matrix& i) : Matrix_base<complex<double>>(r.ndim(), r.mdim()) {
  assert(r.ndim() == i.ndim() && r.mdim() == i.mdim());
  add_real_block(complex<double>(1.0, 0.0), 0, 0, ndim_, mdim_, r.data());
  add_real_block(complex<double>(0.0, 1.0), 0, 0, ndim_, mdim_, i.data());
}


ZMatrix::ZMatrix(const Matrix& r, const complex<double> factor) : Matrix_base<complex<double>>(r.ndim(), r.mdim()) {
  add_real_block(factor, 0, 0, ndim_, mdim_, r.data());
}


#ifdef HAVE_SCALAPACK
ZMatrix::ZMatrix(const DistZMatrix& o) : Matrix_base<complex<double>>(o.ndim(), o.mdim()) {
  setlocal_(o.local());
}
#endif


ZMatrix ZMatrix::operator+(const ZMatrix& o) const {
  ZMatrix out(*this);
  out.ax_plus_y(complex<double>(1.0,0.0), o);
  return out;
}


ZMatrix& ZMatrix::operator+=(const ZMatrix& o) {
  ax_plus_y(complex<double>(1.0,0.0), o);
  return *this;
}


ZMatrix& ZMatrix::operator-=(const ZMatrix& o) {
  ax_plus_y(complex<double>(-1.0,0.0), o);
  return *this;
}


ZMatrix& ZMatrix::operator=(const ZMatrix& o) {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  copy_n(o.data(), ndim_*mdim_, data());
  return *this;
}


ZMatrix& ZMatrix::operator=(ZMatrix&& o) {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  data_ = move(o.data_);
  return *this;
}


ZMatrix ZMatrix::operator-(const ZMatrix& o) const {
  ZMatrix out(*this);
  out.ax_plus_y(complex<double>(-1.0,0.0), o);
  return out;
}


ZMatrix ZMatrix::operator*(const ZMatrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.ndim());
  const int n = o.mdim();
  ZMatrix out(l, n, localized_);

#ifdef HAVE_SCALAPACK
  assert(localized_ == o.localized_);
  if (localized_) {
#endif
    zgemm3m_("N", "N", l, n, m, 1.0, data(), l, o.data(), o.ndim_, 0.0, out.data(), l);
#ifdef HAVE_SCALAPACK
  } else {
    unique_ptr<complex<double>[]> locala = getlocal();
    unique_ptr<complex<double>[]> localb = o.getlocal();
    unique_ptr<complex<double>[]> localc = out.getlocal();
    pzgemm_("N", "N", l, n, m, 1.0, locala.get(), desc_.get(), localb.get(), o.desc_.get(), 0.0, localc.get(), out.desc_.get());
    out.setlocal_(localc);
  }
#endif

  return out;
}


ZMatrix& ZMatrix::operator*=(const ZMatrix& o) {
  *this = *this * o;
  return *this;
}


ZMatrix ZMatrix::operator*(const complex<double>& a) const {
  ZMatrix out(*this);
  out *= a;
  return out;
}


ZMatrix ZMatrix::operator/(const complex<double>& a) const {
  ZMatrix out(*this);
  out /= a;
  return out;
}


ZMatrix& ZMatrix::operator*=(const complex<double>& a) {
  zscal_(ndim_*mdim_, a, data_, 1);
  return *this;
}
ZMatrix& ZMatrix::operator/=(const complex<double>& a) {
  *this *= 1.0/a;
  return *this;
}


ZMatrix ZMatrix::operator%(const ZMatrix& o) const {
  const int l = mdim_;
  const int m = ndim_;
  assert(ndim_ == o.ndim());
  const int n = o.mdim();
  ZMatrix out(l, n, localized_);

#ifdef HAVE_SCALAPACK
  assert(localized_ == o.localized_);
  if (localized_) {
#endif
    zgemm3m_("C", "N", l, n, m, 1.0, data(), m, o.data(), o.ndim_, 0.0, out.data(), l);
#ifdef HAVE_SCALAPACK
  } else {
    unique_ptr<complex<double>[]> locala = getlocal();
    unique_ptr<complex<double>[]> localb = o.getlocal();
    unique_ptr<complex<double>[]> localc = out.getlocal();
    pzgemm_("C", "N", l, n, m, 1.0, locala.get(), desc_.get(), localb.get(), o.desc_.get(), 0.0, localc.get(), out.desc_.get());
    out.setlocal_(localc);
  }
#endif

  return out;
}


ZMatrix ZMatrix::operator^(const ZMatrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.mdim());
  const int n = o.ndim();
  ZMatrix out(l, n, localized_);

#ifdef HAVE_SCALAPACK
  assert(localized_ == o.localized_);
  if (localized_) {
#endif
    zgemm3m_("N", "C", l, n, m, 1.0, data(), ndim_, o.data(), o.ndim_, 0.0, out.data(), l);
#ifdef HAVE_SCALAPACK
  } else {
    unique_ptr<complex<double>[]> locala = getlocal();
    unique_ptr<complex<double>[]> localb = o.getlocal();
    unique_ptr<complex<double>[]> localc = out.getlocal();
    pzgemm_("N", "C", l, n, m, 1.0, locala.get(), desc_.get(), localb.get(), o.desc_.get(), 0.0, localc.get(), out.desc_.get());
    out.setlocal_(localc);
  }
#endif

  return out;
}


ZMatrix ZMatrix::operator/(const ZMatrix& o) const {
  ZMatrix out(*this);
  out /= o;
  return out;
}


ZMatrix& ZMatrix::operator/=(const ZMatrix& o) {
  assert(ndim_ == o.ndim_); assert(mdim_ == o.mdim_);
  auto oiter = o.cbegin();
  for (auto& i : *this) {
    i /= *oiter++;
  }
  return *this;
}


void ZMatrix::diagonalize(double* eig) {
  if (ndim_ != mdim_) throw logic_error("illegal call of ZMatrix::diagonalize(complex<double>*)");

  //assert that matrix is hermitian to ensure real eigenvalues
  assert((*this - *(this->transpose_conjg())).norm()/size() < 1e-10);

  const int n = ndim_;
  int info;
#ifdef HAVE_SCALAPACK
  if (localized_) {
#endif
    unique_ptr<complex<double>[]> work(new complex<double>[n*6]);
    unique_ptr<double[]> rwork(new double[3*ndim_]);
    zheev_("V", "L", n, data(), n, eig, work.get(), n*6, rwork.get(), info);
    mpi__->broadcast(data(), n*n, 0);
#ifdef HAVE_SCALAPACK
  } else {
    const int localrow = get<0>(localsize_);
    const int localcol = get<1>(localsize_);

    unique_ptr<complex<double>[]> coeff(new complex<double>[localrow*localcol]);
    unique_ptr<complex<double>[]> local = getlocal();

    // first compute worksize
    complex<double> wsize;
    const int lrwork = 1 + 9*n + 3*localrow*localcol;
    const int liwork = 7*n + 8*mpi__->npcol() + 2;
    unique_ptr<double[]> rwork(new double[lrwork]);
    unique_ptr<int[]> iwork(new int[liwork]);

    pzheevd_("V", "U", n, local.get(), desc_.get(), eig, coeff.get(), desc_.get(), &wsize, -1, rwork.get(), lrwork, iwork.get(), liwork, info);

    const int lwork = round(max(131072.0, wsize.real()*2.0));
    unique_ptr<complex<double>[]> work(new complex<double>[lwork]);
    pzheevd_("V", "U", n, local.get(), desc_.get(), eig, coeff.get(), desc_.get(), work.get(), lwork, rwork.get(), lrwork, iwork.get(), liwork, info);
    if(info) throw runtime_error("diagonalize failed");

    setlocal_(coeff);
  }
#endif

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
  zgesvd_("A", "A", ndim_, mdim_, cblock, ndim_, S.get(), ublock, ndim_, vblock, mdim_, work.get(), lwork, rwork.get(), info);
  if (info != 0) throw runtime_error("zgesvd failed in ZMatrix::svd");
}


// this is a vector B, and solve AC = B, returns C
shared_ptr<ZMatrix> ZMatrix::solve(shared_ptr<const ZMatrix> A, const int n) const {
  ZMatrix As = *A;
  auto out = make_shared<ZMatrix>(*this);
  assert(n <= out->ndim() && n <= A->ndim() && n <= A->mdim());

  unique_ptr<int[]> ipiv(new int[n]);
  int info;
  zgesv_(n, out->mdim(), As.data(), As.ndim(), ipiv.get(), out->data(), out->ndim(), info);
  if (info) throw runtime_error("ZGESV failed");

  return out;
}


shared_ptr<ZMatrix> ZMatrix::exp(const int deg) const {
  auto out = make_shared<ZMatrix>(ndim_, mdim_, localized_);
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
  auto out = make_shared<ZMatrix>(ndim_, mdim_, localized_);
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
  auto out = make_shared<ZMatrix>(mdim_, ndim_, localized_);
  mytranspose_(data_.get(), ndim_, mdim_, out->data());
  return out;
}


shared_ptr<ZMatrix> ZMatrix::transpose_conjg() const {
  auto out = make_shared<ZMatrix>(mdim_, ndim_, localized_);
  mytranspose_conjg_(data_.get(), ndim_, mdim_, out->data());
  return out;
}

void ZMatrix::antisymmetrize() {
  assert(ndim_ == mdim_);

  shared_ptr<ZMatrix> trans = transpose();
  *this -= *trans;
  *this *= 0.5;
}


void ZMatrix::hermite() {
  assert(ndim_ == mdim_);

  shared_ptr<ZMatrix> trans = transpose_conjg();
  *this += *trans;
  *this *= 0.5;
}


void ZMatrix::purify_unitary() {
  assert(ndim_ == mdim_);
  // Schmidt orthogonalization
  for (int i = 0; i != ndim_; ++i) {
    for (int j = 0; j != i; ++j) {
      const complex<double> a = zdotc_(ndim_, &data_[j*ndim_], 1, &data_[i*ndim_], 1);
      zaxpy_(ndim_, -a, &data_[j*ndim_], 1, &data_[i*ndim_], 1);
    }
    const complex<double> b = 1.0/sqrt(zdotc_(ndim_, &data_[i*ndim_], 1, &data_[i*ndim_], 1));
    zscal_(ndim_, b, &data_[i*ndim_], 1);
  }
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
      const complex<double> ele = (element(j,i) - conj(element(i,j))) * 0.5;
      element(j,i) = ele;
      element(i,j) = - conj(ele);
    }
  }
#endif

}


void ZMatrix::purify_idempotent(const ZMatrix& s) {
  *this = *this * s * *this * 3.0 - *this * s * *this * s * *this * 2.0;
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


void ZMatrix::print(const string name, const size_t size) const {
  cout << "++++ " + name + " ++++" << endl;
  for (int i = 0; i != min(size,ndim_); ++i) {
    for (int j = 0; j != min(size,mdim_); ++j) {
      cout << fixed << setw(30) << setprecision(8) << data_[j * ndim_ + i]  << " ";
    }
    cout << endl;
  }
}

void ZMatrix::copy_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const double* data) {
  for (int i = mdim_i, j = 0; i != mdim_i + mdim ; ++i, ++j) {
    for (int k = ndim_i, l = 0; k != ndim_i + ndim ; ++k, ++l) {
      element(k,i) = a * *(data + j*ndim + l);
    }
  }
}


void ZMatrix::copy_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const unique_ptr<double[]> data) {
  copy_real_block(a, ndim_i, mdim_i, ndim, mdim, data.get());
}


void ZMatrix::copy_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const shared_ptr<const Matrix> data) {
  assert(ndim == data->ndim() && mdim == data->mdim());
  copy_real_block(a, ndim_i, mdim_i, ndim, mdim, data->data());
}


void ZMatrix::copy_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const Matrix& data) {
  assert(ndim == data.ndim() && mdim == data.mdim());
  copy_real_block(a, ndim_i, mdim_i, ndim, mdim, data.data());
}


void ZMatrix::add_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const double* data) {
  for (int i = mdim_i, j = 0; i != mdim_i + mdim ; ++i, ++j) {
    for (int k = ndim_i, l = 0; k != ndim_i + ndim ; ++k, ++l) {
      element(k,i) += a * *(data + j*ndim + l);
    }
  }
}


void ZMatrix::add_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const unique_ptr<double[]> data) {
  add_real_block(a, ndim_i, mdim_i, ndim, mdim, data.get());
}


void ZMatrix::add_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const shared_ptr<const Matrix> data) {
  assert(ndim == data->ndim() && mdim == data->mdim());
  add_real_block(a, ndim_i, mdim_i, ndim, mdim, data->data());
}


void ZMatrix::add_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const Matrix& data) {
  assert(ndim == data.ndim() && mdim == data.mdim());
  add_real_block(a, ndim_i, mdim_i, ndim, mdim, data.data());
}


shared_ptr<Matrix> ZMatrix::get_real_part() const {
  auto out = make_shared<Matrix>(ndim_, mdim_, localized_);
  auto i = cbegin();
  for (auto& o : *out)
    o = real(*i++);
  return out;
}


shared_ptr<Matrix> ZMatrix::get_imag_part() const {
  auto out = make_shared<Matrix>(ndim_, mdim_, localized_);
  auto i = cbegin();
  for (auto& o : *out)
    o = imag(*i++);
  return out;
}


shared_ptr<ZMatrix> ZMatrix::get_conjg() const {
  auto out = make_shared<ZMatrix>(ndim_, mdim_, localized_);
  auto i = cbegin();
  for (auto& o : *out)
    o = conj(*i++);
  return out;
}


#ifdef HAVE_SCALAPACK
shared_ptr<DistZMatrix> ZMatrix::distmatrix() const {
  return make_shared<DistZMatrix>(*this);
}
#else
shared_ptr<const ZMatrix> ZMatrix::distmatrix() const {
  shared_ptr<const ZMatrix> out = shared_from_this();
  return out;
}
#endif


#ifndef HAVE_SCALAPACK
shared_ptr<const ZMatrix> ZMatrix::form_density_rhf(const int n, const int offset) const {
  shared_ptr<const ZMatrix> tmp = this->slice(offset, offset+n);
  auto out = make_shared<const ZMatrix>(*tmp ^ *tmp);
  return out;
}
#endif

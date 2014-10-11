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
#include <src/math/matop.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(ZMatrix)

ZMatrix::ZMatrix(const int n, const int m, const bool loc) : Matrix_base<complex<double>>(n, m, loc) {
}


ZMatrix::ZMatrix(const ZMatrix& o) : Matrix_base<complex<double>>(o) {
}


ZMatrix::ZMatrix(const ZMatView& o) : Matrix_base<complex<double>>(o) {
}


ZMatrix::ZMatrix(ZMatrix&& o) : Matrix_base<complex<double>>(move(o)) {
}


ZMatrix::ZMatrix(const Matrix& r, const Matrix& i) : Matrix_base<complex<double>>(r.ndim(), r.mdim()) {
  assert(r.ndim() == i.ndim() && r.mdim() == i.mdim());
  add_real_block(complex<double>(1.0, 0.0), 0, 0, ndim(), mdim(), r);
  add_real_block(complex<double>(0.0, 1.0), 0, 0, ndim(), mdim(), i);
}


ZMatrix::ZMatrix(const Matrix& r, const complex<double> factor) : Matrix_base<complex<double>>(r.ndim(), r.mdim(), r.localized()) {
  add_real_block(factor, 0, 0, ndim(), mdim(), r);
}


#ifdef HAVE_SCALAPACK
ZMatrix::ZMatrix(const DistZMatrix& o) : Matrix_base<complex<double>>(o.ndim(), o.mdim()) {
  setlocal_(o.local());
}
#endif


ZMatrix ZMatrix::operator/(const ZMatrix& o) const {
  ZMatrix out(*this);
  out /= o;
  return out;
}


ZMatrix& ZMatrix::operator/=(const ZMatrix& o) {
  assert(ndim() == o.ndim()); assert(mdim() == o.mdim());
  auto oiter = o.cbegin();
  for (auto& i : *this) {
    i /= *oiter++;
  }
  return *this;
}


void ZMatrix::diagonalize(VecView eig) {
  if (ndim() != mdim()) throw logic_error("illegal call of ZMatrix::diagonalize(complex<double>*)");
  assert(eig.size() >= ndim());

  // assert that matrix is hermitian to ensure real eigenvalues
  assert(is_hermitian(1.0e-10));

  const int n = ndim();
  int info;
#ifdef HAVE_SCALAPACK
  if (localized_  || n <= blocksize__) {
#endif
    unique_ptr<complex<double>[]> work(new complex<double>[n*6]);
    unique_ptr<double[]> rwork(new double[3*ndim()]);
    zheev_("V", "L", n, data(), n, eig.data(), work.get(), n*6, rwork.get(), info);
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

    pzheevd_("V", "U", n, local.get(), desc_.data(), eig.data(), coeff.get(), desc_.data(), &wsize, -1, rwork.get(), lrwork, iwork.get(), liwork, info);

    const int lwork = round(max(131072.0, wsize.real()*2.0));
    unique_ptr<complex<double>[]> work(new complex<double>[lwork]);
    pzheevd_("V", "U", n, local.get(), desc_.data(), eig.data(), coeff.get(), desc_.data(), work.get(), lwork, rwork.get(), lrwork, iwork.get(), liwork, info);

    setlocal_(coeff);
  }
#endif

  if(info) throw runtime_error("diagonalize failed");

#ifndef NDEBUG
  double max_eig = 0.0;
  double min_eig = 100.0;
  for (int i=0; i!=ndim(); ++i) {
    if (std::abs(eig(i)) > max_eig) max_eig = std::abs(eig(i));
    if (std::abs(eig(i)) < min_eig) min_eig = std::abs(eig(i));
  }
  const double condition_number = max_eig / min_eig;
  if (condition_number > 1.0e8)
    cout << "    - condition number = " << setw(14) << scientific << setprecision(4) << condition_number << " - watch for numerical error in diagonalization" << endl;
#endif
}


tuple<shared_ptr<ZMatrix>, shared_ptr<ZMatrix>> ZMatrix::svd(double* sing) {
  auto U = make_shared<ZMatrix>(ndim(), ndim());
  auto V = make_shared<ZMatrix>(mdim(), mdim());
  const int lwork = 10*max(ndim(), mdim());
  unique_ptr<double[]> rwork(new double[5*max(ndim(), mdim())]);
  unique_ptr<complex<double>[]> work(new complex<double>[lwork]);

  unique_ptr<double[]> S;
  if (!sing) {
    S = unique_ptr<double[]>(new double[min(ndim(), mdim())]);
    sing = S.get();
  }
/*
  SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 $                   WORK, LWORK, RWORK, INFO )
 */
  complex<double>* cblock = data();
  complex<double>* ublock = U->data();
  complex<double>* vblock = V->data();
  int info = 0;
  zgesvd_("A", "A", ndim(), mdim(), cblock, ndim(), sing, ublock, ndim(), vblock, mdim(), work.get(), lwork, rwork.get(), info);
  if (info != 0) throw runtime_error("zgesvd failed in ZMatrix::svd");

  return make_tuple(U, V);
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
  auto out = make_shared<ZMatrix>(ndim(), mdim(), localized_);
  ZMatrix buf(*this);
  assert(ndim() == mdim());

  for (int i = deg; i != 1; --i) {
    const complex<double> inv = 1.0/static_cast<complex<double>>(i);
    buf *= inv;
    for (int j = 0; j != ndim(); ++j) buf(j,j) += 1.0;
    *out = (*this)*buf;
    if (i != 1) buf = *out;
  }
  for (int j = 0; j != ndim(); ++j) out->element(j,j) += 1.0;
  return out;
}


shared_ptr<ZMatrix> ZMatrix::log(const int deg) const {
  auto out = make_shared<ZMatrix>(ndim(), mdim(), localized_);
  ZMatrix buf(*this);
  for (int j = 0; j != ndim(); ++j) buf(j,j) -= 1.0;
  assert(ndim() == mdim());

  for (int i = deg; i != 1; --i) {
    const complex<double> inv = -static_cast<complex<double>>(i-1)/static_cast<complex<double>>(i);
    buf *= inv;
    for (int j = 0; j != ndim(); ++j) buf(j,j) += 1.0;
    *out = (*this)*buf - buf;
    if (i != 1) buf = *out;
  }
  return out;
}


unique_ptr<complex<double>[]> ZMatrix::diag() const {
  if (ndim() != mdim()) throw logic_error("illegal call of ZMatrix::diag()");
  unique_ptr<complex<double>[]> out(new complex<double>[ndim()]);
  for (int i = 0; i != ndim(); ++i) {
    out[i] = element(i,i);
  }
  return move(out);
}


shared_ptr<ZMatrix> ZMatrix::transpose(const complex<double> factor) const {
  auto out = make_shared<ZMatrix>(mdim(), ndim(), localized_);
  blas::transpose(data(), ndim(), mdim(), out->data(), factor);
  return out;
}


shared_ptr<ZMatrix> ZMatrix::transpose_conjg(const complex<double> factor) const {
  auto out = make_shared<ZMatrix>(mdim(), ndim(), localized_);
  blas::transpose_conjg(data(), ndim(), mdim(), out->data(), factor);
  return out;
}


void ZMatrix::purify_unitary() {
  assert(ndim() == mdim());
  // Schmidt orthogonalization
  for (int i = 0; i != ndim(); ++i) {
    for (int j = 0; j != i; ++j) {
      const complex<double> a = blas::dot_product(element_ptr(0,j), ndim(), element_ptr(0,i));
      blas::ax_plus_y_n(-a, element_ptr(0,j), ndim(), element_ptr(0,i));
    }
    const complex<double> b = 1.0/sqrt(blas::dot_product(element_ptr(0,i), ndim(), element_ptr(0,i)));
    for_each(element_ptr(0,i), element_ptr(0,i+1), [&b](complex<double>& a) { a *= b; });
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
  for (int i = 0; i != ndim(); ++i) {
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
  assert(ndim() == mdim());
  const int n = ndim();
#ifndef NDEBUG
  shared_ptr<ZMatrix> ref = this->copy();
#endif
  shared_ptr<ZMatrix> buf = this->clone();
  buf->unit();

  int info;
  unique_ptr<int[]> ipiv(new int[n]);
  zgesv_(n, n, data(), n, ipiv.get(), buf->data(), n, info);
  if (info) throw runtime_error("dsysv failed in ZMatrix::inverse()");

  // check numerical stability of the inversion
#ifndef NDEBUG
  assert((*ref * *buf).is_identity());
#endif

  copy_n(buf->data(), n*n, data());
}


// compute S^{-1/2}
bool ZMatrix::inverse_half(const double thresh) {
  assert(ndim() == mdim());
  const int n = ndim();
#ifndef NDEBUG
  shared_ptr<ZMatrix> ref = this->copy();
#endif
  VectorB vec(n);
  diagonalize(vec);

  for (int i = 0; i != n; ++i) {
    double s = vec(i) > thresh ? 1.0/sqrt(sqrt(vec(i))) : 0.0;
    for_each(element_ptr(0,i), element_ptr(0,i+1), [&s](complex<double>& a) { a *= s; });
  }

#ifndef NDEBUG
  for (int i = 0; i != n; ++i)
    if (vec[i] < thresh) cout << " throwing out " << setprecision(20) << vec[i] << endl;
#endif

  *this = *this ^ *this;

  const bool lindep = std::any_of(vec.begin(), vec.end(), [&thresh] (const double& e) { return e < thresh; });

#ifndef NDEBUG
  // check numerical stability of the inversion - bypassed if we detect linear dependency
  assert((*this % *ref * *this).is_identity() || lindep);
#endif

  return !lindep;
}


shared_ptr<ZMatrix> ZMatrix::tildex(const double thresh) const {
  shared_ptr<ZMatrix> out = this->copy();
  bool nolindep = out->inverse_half(thresh);
  if (!nolindep) {
    // use canonical orthogonalization. Start over
    cout << "    * Using canonical orthogonalization due to linear dependency" << endl << endl;
    out = this->copy();
    VectorB eig(ndim());
    out->diagonalize(eig);
    int m = 0;
    for (int i = 0; i != mdim(); ++i) {
      if (eig[i] > thresh) {
        const double e = 1.0/std::sqrt(eig(i));
        transform(out->element_ptr(0,i), out->element_ptr(0,i+1), out->element_ptr(0,m++), [&e](complex<double> a) { return a*e; });
      }
    }
    out = out->slice_copy(0,m);
  }

  // check numerical stability of the orthogonalization
  assert((*out % *this * *out).is_identity());

  return out;
}


void ZMatrix::copy_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const MatView data) {
  assert(ndim == data.ndim() && mdim == data.mdim());
  for (int i = mdim_i, j = 0; i != mdim_i + mdim ; ++i, ++j)
    transform(&data(0,j), &data(0,j+1), element_ptr(ndim_i,i), [&a](const double& p) { return a*p; });
}


void ZMatrix::add_real_block(const complex<double> a, const int ndim_i, const int mdim_i, const int ndim, const int mdim, const MatView data) {
  assert(ndim == data.ndim() && mdim == data.mdim());
  for (int i = mdim_i, j = 0; i != mdim_i + mdim ; ++i, ++j)
    blas::ax_plus_y_n(a, data.element_ptr(0, j), ndim, element_ptr(ndim_i, i));
}


shared_ptr<Matrix> ZMatrix::get_real_part() const {
  auto out = make_shared<Matrix>(ndim(), mdim(), localized_);
  auto i = cbegin();
  for (auto& o : *out)
    o = real(*i++);
  return out;
}


shared_ptr<Matrix> ZMatrix::get_imag_part() const {
  auto out = make_shared<Matrix>(ndim(), mdim(), localized_);
  auto i = cbegin();
  for (auto& o : *out)
    o = imag(*i++);
  return out;
}


shared_ptr<ZMatrix> ZMatrix::get_conjg() const {
  auto out = make_shared<ZMatrix>(ndim(), mdim(), localized_);
  auto i = cbegin();
  for (auto& o : *out)
    o = conj(*i++);
  return out;
}


bool ZMatrix::is_symmetric(const double thresh) const {
  shared_ptr<ZMatrix> A = copy();
  *A -= *A->transpose();
  const double err = A->norm()/A->size();
#ifndef NDEBUG
  if (100.0*err > thresh)
    cout << scientific << setprecision(2) << "    - ZMatrix symmetry not fully satisfied: error norm/size = " << err << std::endl;
#endif
  return (err < thresh);
}


bool ZMatrix::is_antisymmetric(const double thresh) const {
  shared_ptr<ZMatrix> A = copy();
  *A += *A->transpose();
  const double err = A->norm()/A->size();
#ifndef NDEBUG
  if (100.0*err > thresh)
    cout << scientific << setprecision(2) << "    - ZMatrix antisymmetry not fully satisfied: error norm/size = " << err << std::endl;
#endif
  return (err < thresh);
}


bool ZMatrix::is_hermitian(const double thresh) const {
  shared_ptr<ZMatrix> A = copy();
  *A -= *A->transpose_conjg();
  const double err = A->norm()/A->size();
#ifndef NDEBUG
  if (100.0*err > thresh)
    cout << scientific << setprecision(2) << "    - Hermitian symmetry not fully satisfied: error norm/size = " << err << std::endl;
#endif
  return (err < thresh);
}


bool ZMatrix::is_identity(const double thresh) const {
  shared_ptr<ZMatrix> A = copy();
  shared_ptr<ZMatrix> B = A->clone();
  B->unit();
  *A -= *B;
  const double err = A->norm()/A->size();
#ifndef NDEBUG
  if (100.0*err > thresh)
    cout << scientific << setprecision(2) << "    - Inversion not perfectly accurate: error norm/size = " << err << std::endl;
#endif
  return (err < thresh);
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


shared_ptr<const ZMatrix> ZMatrix::form_density_rhf(const int n, const int offset, const complex<double> scale) const {
  const ZMatView tmp = this->slice(offset, offset+n);
  auto out = make_shared<ZMatrix>(tmp ^ tmp);
  *out *= scale;
  return out;
}

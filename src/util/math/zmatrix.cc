//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zmatrix.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include <iostream>
#include <iomanip>
#include <cmath>
#include <src/util/taskqueue.h>
#include <src/util/constants.h>
#include <src/util/math/zmatrix.h>
#include <src/util/math/matop.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(ZMatrix)

ZMatrix::ZMatrix(const int n, const int m, const bool loc) : Matrix_base<complex<double>>(n, m, loc), std::enable_shared_from_this<bagel::ZMatrix>() {
}


ZMatrix::ZMatrix(const ZMatrix& o) : Matrix_base<complex<double>>(o), std::enable_shared_from_this<bagel::ZMatrix>() {
}


ZMatrix::ZMatrix(const ZMatView& o) : Matrix_base<complex<double>>(o), std::enable_shared_from_this<bagel::ZMatrix>() {
}


ZMatrix::ZMatrix(ZMatrix&& o) : Matrix_base<complex<double>>(move(o)), std::enable_shared_from_this<bagel::ZMatrix>() {
}


ZMatrix::ZMatrix(const Matrix& r, const Matrix& i) : Matrix_base<complex<double>>(r.ndim(), r.mdim()), std::enable_shared_from_this<bagel::ZMatrix>() {
  assert(r.ndim() == i.ndim() && r.mdim() == i.mdim());
  add_real_block(complex<double>(1.0, 0.0), 0, 0, ndim(), mdim(), r);
  add_real_block(complex<double>(0.0, 1.0), 0, 0, ndim(), mdim(), i);
}


ZMatrix::ZMatrix(const Matrix& r, const complex<double> factor) : Matrix_base<complex<double>>(r.ndim(), r.mdim(), r.localized()), std::enable_shared_from_this<bagel::ZMatrix>() {
  add_real_block(factor, 0, 0, ndim(), mdim(), r);
}


#ifdef HAVE_SCALAPACK
ZMatrix::ZMatrix(const DistZMatrix& o) : Matrix_base<complex<double>>(o.ndim(), o.mdim()), std::enable_shared_from_this<bagel::ZMatrix>() {
  setlocal_(o.local());
}
#endif


ZMatView ZMatrix::slice(const int mstart, const int mend) {
  assert(mstart >= 0 && mend <= mdim());
  auto low = {0, mstart};
  auto up  = {static_cast<int>(ndim()), mend};
  return ZMatView(btas::make_rwview(this->range().slice(low, up), this->storage()), localized_);
}


const ZMatView ZMatrix::slice(const int mstart, const int mend) const {
  assert(mstart >= 0 && mend <= mdim());
  auto low = {0, mstart};
  auto up  = {static_cast<int>(ndim()), mend};
  return ZMatView(btas::make_rwview(this->range().slice(low, up), this->storage()), localized_);
}


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
//assert(is_hermitian(1.0e-10));

  const int n = ndim();
  int info;
#ifdef HAVE_SCALAPACK
  if (localized_  || n <= 10*blocksize__) {
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
    const complex<double> b = 1.0/std::sqrt(blas::dot_product(element_ptr(0,i), ndim(), element_ptr(0,i)));
    for_each(element_ptr(0,i), element_ptr(0,i+1), [&b](complex<double>& a) { a *= b; });
  }
}


// in-place matrix inverse (practically we use buffer area)
void ZMatrix::inverse() {
  assert(ndim() == mdim());
  const int n = ndim();
  shared_ptr<ZMatrix> buf = this->clone();
  buf->unit();

  int info;
  unique_ptr<int[]> ipiv(new int[n]);
  zgesv_(n, n, data(), n, ipiv.get(), buf->data(), n, info);
  if (info) throw runtime_error("dsysv failed in ZMatrix::inverse()");

  copy_n(buf->data(), n*n, data());
}


// compute S^{-1/2}
bool ZMatrix::inverse_half(const double thresh) {
  assert(ndim() == mdim());
  const int n = ndim();
  VectorB vec(n);
  diagonalize(vec);

  for (int i = 0; i != n; ++i) {
    double s = vec(i) > thresh ? 1.0/std::sqrt(std::sqrt(vec(i))) : 0.0;
    for_each(element_ptr(0,i), element_ptr(0,i+1), [&s](complex<double>& a) { a *= s; });
  }
  *this = *this ^ *this;
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
void ZMatrix::sqrt() {
  assert(ndim() == mdim());
  const int n = ndim();
  VectorB vec(n);
#ifdef HAVE_SCALAPACK
  if (localized_) {
#endif
    diagonalize(vec);

    for (int i = 0; i != n; ++i) {
      if (vec(i) < -numerical_zero__) throw runtime_error("Matrix::sqrt() called, but this matrix is not positive definite");
      blas::scale_n(std::sqrt(std::sqrt(fabs(vec(i)))), element_ptr(0,i), n);
    }

    *this = *this ^ *this;
#ifdef HAVE_SCALAPACK
  }
  else {
    unique_ptr<double[]> scal(new double[n]);
    shared_ptr<DistZMatrix> dist = distmatrix();
    dist->diagonalize(vec);
    for (int i = 0; i != n; ++i)
      scal[i] = std::sqrt(std::sqrt(vec(i)));
    dist->scale(scal.get());
    *this = *(*dist ^ *dist).matrix();
  }
#endif
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
  shared_ptr<ZMatrix> tmp = transpose();
  *tmp -= *this;
  return tmp->rms() < thresh;
}


bool ZMatrix::is_antisymmetric(const double thresh) const {
  shared_ptr<ZMatrix> tmp = transpose();
  *tmp += *this;
  return tmp->rms() < thresh;
}


bool ZMatrix::is_hermitian(const double thresh) const {
  shared_ptr<ZMatrix> tmp = transpose_conjg();
  *tmp -= *this;
  return tmp->rms() < thresh;
}


bool ZMatrix::is_identity(const double thresh) const {
  shared_ptr<ZMatrix> tmp = clone();
  tmp->unit();
  *tmp -= *this;
  return tmp->rms() < thresh;
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


void ZMatrix::print (const string tag, int len) const {
  if (tag != "")
    cout << endl << "  ++ " << tag << " ++" << endl << endl;

  int len_n, len_m;

  if (len == 0 || (len > mdim() && len > ndim())) {
    len_n = ndim();
    len_m = mdim();
  } else {
    len_n = len_m = len;
  }

  for (int i = 0; i != len_m / 6; ++i) {
    cout << setw(6) << " ";
    for (int k = 0; k != 6; ++k)
      cout << setw(30) << i * 6 + k;
    cout << endl;
    for (int j = 0; j != len_n; ++j) {
      cout << setw(6) << j;
      for (int k = 0; k != 6; ++k)
        cout << setw(30) << setprecision(10) << *element_ptr(j, i * 6 + k);
      cout << endl;
    }
    cout << endl;
  }

  if (len_m % 6) {
    int i = len_m / 6;
    cout << setw(6) << " ";
    for (int k = 0; k != len_m % 6; ++k)
      cout << setw(30) << i * 6 + k;
    cout << endl;

    for (int j = 0; j != len_n; ++j) {
      cout << setw(6) << j;
      for (int k = 0; k != len_m % 6; ++k)
        cout << setw(30) << setprecision(10) << *element_ptr(j, i * 6 + k);
      cout << endl;
    }
    cout << endl;
  }
}



shared_ptr<const ZMatrix> ZMatrix::form_density_rhf(const int n, const int offset, const complex<double> scale) const {
  const ZMatView tmp = this->slice(offset, offset+n);
  auto out = make_shared<ZMatrix>(tmp ^ tmp);
  *out *= scale;
  return out;
}

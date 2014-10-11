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
#include <src/math/matop.h>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <src/util/taskqueue.h>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(Matrix)

Matrix::Matrix(const int n, const int m , const bool loc) : Matrix_base<double>(n,m,loc) {
}


Matrix::Matrix(const Matrix& o) : Matrix_base<double>(o) {
}

Matrix::Matrix(const MatView& o) : Matrix_base<double>(o) {
}

Matrix::Matrix(Matrix&& o) : Matrix_base<double>(move(o)) {
}


#ifdef HAVE_SCALAPACK
Matrix::Matrix(const DistMatrix& o) : Matrix_base<double>(o.ndim(), o.mdim()) {
  setlocal_(o.local());
}
#endif


Matrix Matrix::operator/(const Matrix& o) const {
  Matrix out(*this);
  out /= o;
  return out;
}


Matrix& Matrix::operator/=(const Matrix& o) {
  assert(ndim() == o.ndim()); assert(mdim() == o.mdim());
  auto oiter = o.cbegin();
  for (auto& i : *this) {
    i /= *oiter++;
  }
  return *this;
}


void Matrix::diagonalize(VecView eig) {
  if (ndim() != mdim()) throw logic_error("illegal call of Matrix::diagonalize(double*)");
  assert(eig.size() >= ndim());
  assert(test_symmetric(1.0e-10));
  const int n = ndim();
  int info;

  // assume that the matrix is symmetric
  // the leading order (nbasis supplied)
#ifdef HAVE_SCALAPACK
  if (localized_ || n <= blocksize__) {
#endif
    unique_ptr<double[]> work(new double[n*6]);
    dsyev_("V", "L", n, data(), n, eig.data(), work.get(), n*6, info);
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
    pdsyevd_("V", "U", n, local.get(), desc_.data(), eig.data(), coeff.get(), desc_.data(), &wsize, -1, &liwork, 1, info);
    unique_ptr<int[]> iwork(new int[liwork]);
    wsize =  max(131072.0, wsize*2.0);

    const int lwork = round(wsize);
    unique_ptr<double[]> work(new double[lwork]);
    pdsyevd_("V", "U", n, local.get(), desc_.data(), eig.data(), coeff.get(), desc_.data(), work.get(), lwork, iwork.get(), liwork, info);
    setlocal_(coeff);
  }
#endif

  if (info) throw runtime_error("dsyev/pdsyevd failed in Matrix");
}


tuple<shared_ptr<Matrix>, shared_ptr<Matrix>> Matrix::svd(double* sing) {
  auto U = make_shared<Matrix>(ndim(), ndim());
  auto V = make_shared<Matrix>(mdim(), mdim());

  const int lwork = 10*max(ndim(), mdim());
  unique_ptr<double[]> work(new double[lwork]);

  // If singular values are not requested
  unique_ptr<double[]> S;
  if (!sing) {
    S = unique_ptr<double[]>(new double[min(ndim(), mdim())]);
    sing = S.get();
  }
/*
  SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
 $                   WORK, LWORK, INFO )
 */
  double* cblock = data();
  double* ublock = U->data();
  double* vblock = V->data();
  int info = 0;
  dgesvd_("A", "A", ndim(), mdim(), cblock, ndim(), sing, ublock, ndim(), vblock, mdim(), work.get(), lwork, info);
  if (info != 0) throw runtime_error("dgesvd failed in Matrix::svd");

  return make_tuple(U,V);
}


shared_ptr<Matrix> Matrix::exp(const int deg) const {
  auto out = make_shared<Matrix>(ndim(), mdim());
  Matrix buf(*this);
  assert(ndim() == mdim());

  for (int i = deg; i != 1; --i) {
    const double inv = 1.0/static_cast<double>(i);
    buf *= inv;
    for (int j = 0; j != ndim(); ++j) buf(j,j) += 1.0;
    *out = (*this)*buf;
    if (i != 1) buf = *out;
  }
  for (int j = 0; j != ndim(); ++j) out->element(j,j) += 1.0;
  return out;
}


shared_ptr<Matrix> Matrix::log(const int deg) const {
  auto out = make_shared<Matrix>(ndim(), mdim());
  Matrix buf(*this);
  for (int j = 0; j != ndim(); ++j) buf(j,j) -= 1.0;
  assert(ndim() == mdim());

  for (int i = deg; i != 1; --i) {
    const double inv = -static_cast<double>(i-1)/static_cast<double>(i);
    buf *= inv;
    for (int j = 0; j != ndim(); ++j) buf(j,j) += 1.0;
    *out = (*this)*buf - buf;
    if (i != 1) buf = *out;
  }
  return out;
}


shared_ptr<Matrix> Matrix::transpose(const double factor) const {
  auto out = make_shared<Matrix>(mdim(), ndim(), localized_);
  blas::transpose(data(), ndim(), mdim(), out->data(), factor);
  return out;
}


void Matrix::purify_unitary() {
  assert(ndim() == mdim());
  for (int i = 0; i != ndim(); ++i) {
    for (int j = 0; j != i; ++j) {
      const double a = blas::dot_product(element_ptr(0,i), ndim(), element_ptr(0,j));
      blas::ax_plus_y_n(-a, element_ptr(0,j), ndim(), element_ptr(0,i));
    }
    const double b = 1.0/std::sqrt(blas::dot_product(element_ptr(0,i), ndim(), element_ptr(0,i)));
    for_each(element_ptr(0,i), element_ptr(0,i+1), [&b](double& a) { a *= b; });
  }
}


void Matrix::purify_redrotation(const int nclosed, const int nact, const int nvirt) {
  assert(ndim() == mdim() && nclosed + nact + nvirt == ndim());
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
      const double ele = (element(j,i) - element(i,j)) * 0.5;
      element(j,i) = ele;
      element(i,j) = -ele;
    }
  }
}


void Matrix::purify_idempotent(const Matrix& s) {
  *this = *this * s * *this * 3.0 - *this * s * *this * s * *this * 2.0;
}


// in-place matrix inverse (practically we use buffer area)
void Matrix::inverse() {
  assert(ndim() == mdim());
  const int n = ndim();
#ifndef NDEBUG
  shared_ptr<Matrix> ref = this->copy();
#endif
  shared_ptr<Matrix> buf = this->clone();
  buf->unit();

  int info;
  unique_ptr<int[]> ipiv(new int[n]);
  dgesv_(n, n, data(), n, ipiv.get(), buf->data(), n, info);
  if (info) throw runtime_error("dgesv failed in Matrix::inverse()");

  // check numerical stability of the inversion
  assert((*ref * *buf).test_unit());

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
  assert(ndim() == mdim());
#ifndef NDEBUG
  shared_ptr<Matrix> ref = this->copy();
#endif
  const int n = ndim();
  VectorB vec(n);
  diagonalize(vec);

  for (int i = 0; i != n; ++i) {
    double s = vec(i) > thresh ? 1.0/std::sqrt(vec[i]) : 0.0;
    for_each(element_ptr(0,i), element_ptr(0,i+1), [&s](double& a){ a*= s; });
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

  // check numerical stability of the inversion
  assert((*this * *ref).test_unit());
#endif

  return rm.empty();
}


// compute S^{-1/2}
bool Matrix::inverse_half(const double thresh) {
  assert(ndim() == mdim());
  const int n = ndim();
  VectorB vec(n);
#ifndef NDEBUG
  shared_ptr<Matrix> ref = this->copy();
#endif

#ifdef HAVE_SCALAPACK
  if (localized_) {
#endif
    diagonalize(vec);
    for (int i = 0; i != n; ++i) {
      double s = vec(i) > thresh ? 1.0/std::sqrt(std::sqrt(vec(i))) : 0.0;
      for_each(element_ptr(0,i), element_ptr(0,i+1), [&s](double& a){ a*= s; });
    }
    *this = *this ^ *this;
#ifdef HAVE_SCALAPACK
  } else {
    unique_ptr<double[]> scal(new double[n]);
    shared_ptr<DistMatrix> dist = distmatrix();
    dist->diagonalize(vec);
    for (int i = 0; i != n; ++i)
      scal[i] = vec(i) > thresh ? 1.0/std::sqrt(std::sqrt(vec(i))) : 0.0;
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

#ifndef NDEBUG
  // check numerical stability of the inversion - bypassed if we detect linear dependency
  assert((*this % *ref * *this).test_unit() || !rm.empty());
#endif

  return rm.empty();
}

shared_ptr<Matrix> Matrix::tildex(const double thresh) const {
  shared_ptr<Matrix> out = this->copy();
  bool nolindep = out->inverse_half(thresh);
  if (!nolindep) {
    // use canonical orthogonalization. Start over
    cout << "    * Using canonical orthogonalization due to linear dependency" << endl << endl;
    out = this->copy();
    VectorB eig(ndim());
    out->diagonalize(eig);
    int m = 0;
    for (int i = 0; i != mdim(); ++i) {
      if (eig(i) > thresh) {
        const double e = 1.0/std::sqrt(eig(i));
        transform(out->element_ptr(0,i), out->element_ptr(0,i+1), out->element_ptr(0,m++), [&e](double a){ return a*e; });
      }
    }
    out = out->slice_copy(0,m);
  }

  // check numerical stability of the orthogonalization
  assert((*out % *this * *out).test_unit());

  return out;
}

// compute Hermitian square root, S^{1/2}
void Matrix::sqrt() {
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
    shared_ptr<DistMatrix> dist = distmatrix();
    dist->diagonalize(vec);
    for (int i = 0; i != n; ++i)
      scal[i] = std::sqrt(std::sqrt(vec(i)));
    dist->scale(scal.get());
    *this = *(*dist ^ *dist).matrix();
  }
#endif
}

// CAUTION: assumes no orbital is rotated twice
void Matrix::rotate(vector<tuple<int, int, double>>& rotations) {
  if (rotations.size() > 6*resources__->max_num_threads()) {
    TaskQueue<function<void(void)>> tq(rotations.size());
    for (auto& irot : rotations)
      tq.emplace_back( [this,irot] { rotate(get<0>(irot), get<1>(irot), get<2>(irot)); } );

    tq.compute();
  }
  else {
    for_each(rotations.begin(), rotations.end(), [this] (tuple<int, int, double> irot) { rotate(get<0>(irot), get<1>(irot), get<2>(irot)); });
  }
}


bool Matrix::test_symmetric(const double thresh) const {
  shared_ptr<Matrix> A = copy();
  *A -= *A->transpose();
  return (A->rms() < thresh);
}


bool Matrix::test_antisymmetric(const double thresh) const {
  shared_ptr<Matrix> A = copy();
  *A += *A->transpose();
  return (A->rms() < thresh);
}


bool Matrix::test_unit(const double thresh) const {
  shared_ptr<Matrix> A = copy();
  shared_ptr<Matrix> B = A->clone();
  B->unit();
  *A -= *B;
  return (A->rms() < thresh);
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


shared_ptr<const Matrix> Matrix::form_density_rhf(const int n, const int offset) const {
  const MatView tmp = this->slice(offset, offset+n);
  auto out = make_shared<Matrix>(tmp ^ tmp);
  *out *= 2.0;
  return out;
}


shared_ptr<Matrix> Matrix::swap_columns(const int i, const int iblock, const int j, const int jblock) const {
    assert(j > i);
    assert(i + iblock - 1 < j);
    assert(j + jblock <= mdim());
    cout << " swap columns " << i << " through " << i + iblock << " with " << j << " through " << j + jblock << endl;
    auto out = this->clone();
    out->copy_block(0, 0, ndim(), i, this->slice(0, i));
    out->copy_block(0, i, ndim(), jblock, this->slice(j, j + jblock));
    out->copy_block(0, i + jblock, ndim(), j - (i+iblock), this->slice(i+iblock, j));
    out->copy_block(0, j + jblock - iblock, ndim(), iblock, this->slice(i, i+iblock));
    out->copy_block(0, j + jblock, ndim(), mdim()-(j+jblock), this->slice(j+jblock, mdim()));
    return out;
}

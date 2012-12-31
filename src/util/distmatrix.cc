//
// BAGEL - Parallel electron correlation program.
// Filename: distmatrix.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <stdexcept>
#include <src/util/matrix.h>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;


#ifdef HAVE_SCALAPACK

DistMatrix::DistMatrix(const int n, const int m) : ndim_(n), mdim_(m), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(mpi__->numroc(ndim_, mdim_)) {
  local_ = unique_ptr<double[]>(new double[size()]);
  zero();
}


DistMatrix::DistMatrix(const DistMatrix& o) : ndim_(o.ndim_), mdim_(o.mdim_), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(mpi__->numroc(ndim_, mdim_)) {
  local_ = unique_ptr<double[]>(new double[size()]);
  *this = o;
}


DistMatrix::DistMatrix(const Matrix& o) : ndim_(o.ndim()), mdim_(o.mdim()), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(mpi__->numroc(ndim_, mdim_)) {
  local_ = o.getlocal();
}


DistMatrix DistMatrix::operator*(const DistMatrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.ndim_);
  const int n = o.mdim_;

  DistMatrix out(l, n);
  pdgemm_("N", "N", l, n, m, 1.0, local_.get(), desc_.get(), o.local_.get(), o.desc_.get(), 0.0, out.local_.get(), out.desc_.get());
  return out;
}


DistMatrix& DistMatrix::operator*=(const DistMatrix& o) {
  *this = *this * o;
  return *this;
}


DistMatrix DistMatrix::operator%(const DistMatrix& o) const {
  const int l = mdim_;
  const int m = ndim_;
  assert(ndim_ == o.ndim_);
  const int n = o.mdim_;

  DistMatrix out(l, n);
  pdgemm_("T", "N", l, n, m, 1.0, local_.get(), desc_.get(), o.local_.get(), o.desc_.get(), 0.0, out.local_.get(), out.desc_.get());
  return out;
}


DistMatrix DistMatrix::operator^(const DistMatrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.mdim_);
  const int n = o.ndim_;

  DistMatrix out(l, n);
  pdgemm_("N", "T", l, n, m, 1.0, local_.get(), desc_.get(), o.local_.get(), o.desc_.get(), 0.0, out.local_.get(), out.desc_.get());
  return out;
}


void DistMatrix::diagonalize(double* eig) {
  if (ndim_ != mdim_) throw logic_error("illegal call of DistMatrix::diagonalize(double*)"); 
  const int n = ndim_;

  DistMatrix tmp(*this);

  // first compute worksize
  double wsize;
  int liwork = 1;
  int info;
  pdsyevd_("V", "U", n, local_.get(), desc_.get(), eig, tmp.local_.get(), tmp.desc_.get(), &wsize, -1, &liwork, 1, info);
  unique_ptr<int[]> iwork(new int[liwork]);
  wsize =  max(131072.0, wsize*2.0);

  const int lwork = round(wsize);
  unique_ptr<double[]> work(new double[lwork]);
  pdsyevd_("V", "U", n, local_.get(), desc_.get(), eig, tmp.local_.get(), tmp.desc_.get(), work.get(), lwork, iwork.get(), liwork, info);
  if (info) throw runtime_error("dsyev failed in DistMatrix");

  *this = tmp;
}


double DistMatrix::ddot(const DistMatrix& o) const {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  double sum = size() ? ddot_(size(), local_.get(), 1, o.local_.get(), 1) : 0.0;
  mpi__->allreduce(&sum, 1);
  return sum;
}


shared_ptr<Matrix> DistMatrix::matrix() const {
  shared_ptr<Matrix> out(new Matrix(*this));
  return out;
}

#endif

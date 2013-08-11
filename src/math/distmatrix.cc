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


#include <stdexcept>
#include <src/math/matrix.h>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;


#ifdef HAVE_SCALAPACK

DistMatrix::DistMatrix(const int n, const int m) : DistMatrix_base<double>(n,m) {}


DistMatrix::DistMatrix(const DistMatrix& o) : DistMatrix_base<double>(o) {}


DistMatrix::DistMatrix(const Matrix& o) : DistMatrix_base<double>(o.ndim(), o.mdim()) {
  copy_n(o.getlocal().get(), size(), local_.get());
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
  if (info) throw runtime_error("pdsyevd failed in DistMatrix");

  // seems MKL does not broadcast for tiny matrices..
  if (n <= blocksize__) mpi__->broadcast(eig, n, 0);

  *this = tmp;
}


shared_ptr<Matrix> DistMatrix::matrix() const {
  return make_shared<Matrix>(*this);
}


shared_ptr<const DistMatrix> DistMatrix::form_density_rhf(const int nocc, const int off) const {
  const int l = ndim_;
  const int n = ndim_;
  auto out = make_shared<DistMatrix>(l, n);
  pdgemm_("N", "T", l, n, nocc, 2.0, local_.get(), 1, 1+off, desc_.get(), local_.get(), 1, 1+off, desc_.get(), 0.0, out->local_.get(), 1, 1, out->desc_.get());
  return out;
}

#endif

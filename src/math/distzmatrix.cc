//
// BAGEL - Parallel electron correlation program.
// Filename: distzmatrix.cc
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
#include <src/math/zmatrix.h>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;


#ifdef HAVE_SCALAPACK

DistZMatrix::DistZMatrix(const int n, const int m) : DistMatrix_base<std::complex<double>>(n,m) {}


DistZMatrix::DistZMatrix(const DistZMatrix& o) : DistMatrix_base<std::complex<double>>(o) {}


DistZMatrix::DistZMatrix(const ZMatrix& o) : DistMatrix_base<std::complex<double>>(o.ndim(), o.mdim()) {
  copy_n(o.getlocal().get(), size(), local_.get());
}


DistZMatrix DistZMatrix::operator*(const DistZMatrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.ndim_);
  const int n = o.mdim_;

  DistZMatrix out(l, n);
  pzgemm_("N", "N", l, n, m, 1.0, local_.get(), desc_.get(), o.local_.get(), o.desc_.get(), 0.0, out.local_.get(), out.desc_.get());
  return out;
}


DistZMatrix& DistZMatrix::operator*=(const DistZMatrix& o) {
  *this = *this * o;
  return *this;
}


DistZMatrix DistZMatrix::operator%(const DistZMatrix& o) const {
  const int l = mdim_;
  const int m = ndim_;
  assert(ndim_ == o.ndim_);
  const int n = o.mdim_;

  DistZMatrix out(l, n);
  pzgemm_("C", "N", l, n, m, 1.0, local_.get(), desc_.get(), o.local_.get(), o.desc_.get(), 0.0, out.local_.get(), out.desc_.get());
  return out;
}


DistZMatrix DistZMatrix::operator^(const DistZMatrix& o) const {
  const int l = ndim_;
  const int m = mdim_;
  assert(mdim_ == o.mdim_);
  const int n = o.ndim_;

  DistZMatrix out(l, n);
  pzgemm_("N", "C", l, n, m, 1.0, local_.get(), desc_.get(), o.local_.get(), o.desc_.get(), 0.0, out.local_.get(), out.desc_.get());
  return out;
}


void DistZMatrix::diagonalize(double* eig) {
  if (ndim_ != mdim_) throw logic_error("illegal call of DistZMatrix::diagonalize(double*)");
  const int n = ndim_;
  const int localrow = get<0>(localsize_);
  const int localcol = get<1>(localsize_);

  DistZMatrix tmp(*this);

  // first compute worksize
  int info;
  std::complex<double> wsize;
  const int lrwork = 1 + 9*n + 3*localrow*localcol;
  const int liwork = 7*n + 8*mpi__->npcol() + 2;
  unique_ptr<double[]> rwork(new double[lrwork]);
  unique_ptr<int[]> iwork(new int[liwork]);
  pzheevd_("V", "U", n, local_.get(), desc_.get(), eig, tmp.local_.get(), tmp.desc_.get(), &wsize, -1, rwork.get(), lrwork, iwork.get(), liwork, info);

  const int lwork = round(max(131072.0, wsize.real()*2.0));
  unique_ptr<complex<double>[]> work(new complex<double>[lwork]);
  pzheevd_("V", "U", n, local_.get(), desc_.get(), eig, tmp.local_.get(), tmp.desc_.get(), work.get(), lwork, rwork.get(), lrwork, iwork.get(), liwork, info);
  if (info) throw runtime_error("pzheevd failed in DistZMatrix");

  *this = tmp;
}


complex<double> DistZMatrix::zdotc(const DistZMatrix& o) const {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  complex<double> sum = size() ? zdotc_(size(), local_.get(), 1, o.local_.get(), 1) : std::complex<double>(0.0,0.0);
  mpi__->allreduce(&sum, 1);
  return sum;
}


shared_ptr<ZMatrix> DistZMatrix::matrix() const {
  return make_shared<ZMatrix>(*this);
}


shared_ptr<const DistZMatrix> DistZMatrix::form_density_rhf(const int nocc, const int off) const {
  const int l = ndim_;
  const int n = ndim_;
  auto out = make_shared<DistZMatrix>(l, n);
  pzgemm_("N", "C", l, n, nocc, 1.0, local_.get(), 1, 1+off, desc_.get(), local_.get(), 1, 1+off, desc_.get(), 0.0, out->local_.get(), 1, 1, out->desc_.get());
  return out;
}

#endif

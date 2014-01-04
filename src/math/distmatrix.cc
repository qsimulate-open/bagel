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
#include <src/parallel/accrequest.h>

using namespace std;
using namespace bagel;


#ifdef HAVE_SCALAPACK

DistMatrix::DistMatrix(const int n, const int m) : DistMatrix_base<double>(n,m) {}


DistMatrix::DistMatrix(const DistMatrix& o) : DistMatrix_base<double>(o) {}


DistMatrix::DistMatrix(DistMatrix&& o) : DistMatrix_base<double>(move(o)) {}


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


DistMatrix& DistMatrix::operator=(const DistMatrix& o) {
  assert(size() == o.size());
  copy_n(o.local_.get(), size(), local_.get());
  return *this;
}


DistMatrix& DistMatrix::operator=(DistMatrix&& o) {
  assert(size() == o.size());
  assert(localsize_ == o.localsize_);
  assert(desc_[0] == o.desc_[0]);
  local_ = move(o.local_);
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

struct RotateTask {
  shared_ptr<const Matrix> buf;
  const int srq;
  const int rrq;
  double* const target;
  pair<double, double> factors;

  RotateTask(shared_ptr<const Matrix> b, const int s, const int r, double* const t, pair<double,double> f):
    buf(b), srq(s), rrq(r), target(t), factors(f) {}
};

// Caution: assumes no repetition of rotation indices (each orbital is involved in at most one rotation)
void DistMatrix::rotate(vector<tuple<int, int, double>> rotations) {
  const int localrow = get<0>(localsize_);
  const int localcol = get<1>(localsize_);

  const int mypcol = mpi__->mypcol();
  const int myprow = mpi__->myprow();

  vector<tuple<int, int, double>> local_rotations;
  vector<RotateTask> rot_tasks;

  // Better way to do this?
  const int rotate_tag = 1 << 8;

  for (auto& irot : rotations) {
    int i = get<0>(irot);
    int j = get<1>(irot);
    int ipcol, jpcol, ioffset, joffset;
    tie(ipcol, ioffset) = locate_column(get<0>(irot));
    tie(jpcol, joffset) = locate_column(get<1>(irot));

    const double gamma = get<2>(irot);

    if ( (ipcol == mypcol) != (jpcol == mypcol) ) {
      const bool iloc = (ipcol == mypcol);

      const int localoffset = iloc ? ioffset : joffset;
      double* const localdata = local_.get() + localoffset * localrow;

      const int remoteoffset = ( iloc ? joffset : ioffset ) * localrow;
      const int remoterank = mpi__->pnum( myprow, ( iloc ? jpcol : ipcol ) );

      auto sbuf = make_shared<Matrix>(localrow, 2);
      copy_n(localdata, localrow, sbuf->element_ptr(0, (iloc?0:1)));

      // send/receive
      const int sendtag = rotate_tag + ( iloc ? i + j*ndim() : j + i*ndim() );
      const int recvtag = rotate_tag + ( iloc ? j + i*ndim() : i + j*ndim() );
      const int srq = mpi__->request_send(sbuf->element_ptr(0,(iloc?0:1)), localrow, remoterank, sendtag);
      const int rrq = mpi__->request_recv(sbuf->element_ptr(0,(iloc?1:0)), localrow, remoterank, recvtag);
      pair<double, double> fac = iloc ? make_pair(cos(gamma), sin(gamma)) : make_pair(-sin(gamma), cos(gamma));
      rot_tasks.emplace_back(sbuf, srq, rrq, localdata, fac);
    }
    else if ( (ipcol == mypcol) && (jpcol == mypcol) ) {
      local_rotations.emplace_back(ioffset, joffset, gamma);
    }
  }

  //rotate locally
  for (auto& irot : local_rotations) {
    const int ioffset = get<0>(irot);
    const int joffset = get<1>(irot);
    const double gamma = get<2>(irot);

    double* const idata = local_.get() + ioffset * localrow;
    double* const jdata = local_.get() + joffset * localrow;

    drot_(localrow, idata, 1, jdata, 1, cos(gamma), sin(gamma));
  }

  for (auto& r : rot_tasks) {
    mpi__->wait(r.rrq);
    const double a = r.factors.first;
    const double b = r.factors.second;
    transform(r.buf->element_ptr(0,0), r.buf->element_ptr(0,1), r.buf->element_ptr(0,1), r.target,
      [&a, &b] (const double& p, const double& q) { return a*p + b*q; });
    mpi__->wait(r.srq);
  }
}


void DistMatrix::rotate(const int i, const int j, const double gamma) {
  rotate(vector<tuple<int, int, double>>(1, make_tuple(i, j, gamma)));
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

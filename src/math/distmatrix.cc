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

// Caution: assumes no repetition of rotation indices (each orbital is involved in at most one rotation)
void DistMatrix::rotate(vector<tuple<int, int, double>> rotations) {
  const int localrow = get<0>(localsize_);
  const int localcol = get<1>(localsize_);

  const int nblock = localrow/blocksize__;
  const int mblock = localcol/blocksize__;

  const int nstride = blocksize__ * mpi__->nprow();
  const int mstride = blocksize__ * mpi__->npcol();

  const int myprow = mpi__->myprow();
  const int mypcol = mpi__->mypcol();

  vector<mutex> mt(localcol);
  for (auto &imt : mt) imt.lock();

  auto sender = make_shared<SendRequest>();
  auto accumer = make_shared<AccRequest>(local_.get(), &mt);

  for (auto& irot : rotations) {
    const int i = get<0>(irot);
    const int j = get<1>(irot);
    const double gamma = get<2>(irot);

    int ipcol, ioffset;
    tie(ipcol, ioffset) = locate_column(i);

    int jpcol, joffset;
    tie(jpcol, joffset) = locate_column(j);

    const bool iloc = (ipcol == mypcol);
    const bool jloc = (jpcol == mypcol);

    if ( iloc && jloc ) {
      // rotate locally
      double* const idata = local_.get() + ioffset * localrow;
      double* const jdata = local_.get() + joffset * localrow;

      drot_(localrow, idata, 1, jdata, 1, cos(gamma), sin(gamma));

      mt[ioffset].unlock();
      mt[joffset].unlock();
    }
    else if ( iloc || jloc ) {
      const int localoffset = iloc ? ioffset : joffset;
      double* const localdata = local_.get() + localoffset * localrow;

      const int remoteoffset = ( iloc ? joffset : ioffset ) * localrow;
      const int remoterank = mpi__->pnum( myprow, ( iloc ? ipcol : jpcol ) );

      unique_ptr<double[]> sbuf(new double[localrow]);
      const double sendfactor = ( iloc ? -sin(gamma) : sin(gamma) ); // double check
      daxpy_(localrow, sendfactor, ldata, 1, sbuf.get(), 1);

      sender->request_send(move(sbuf), localrow, remoterank, remoteoffset);

      transform(localdata, localdata + localrow, localdata, [](double a) { return cos(gamma) * a; });

      // ready to receive now
      mt[ioffset].unlock();
      mt[joffset].unlock();
    }

#ifndef USE_SERVER_THREAD
    sender->flush();
    accumer->flush();
#endif
  }

  // terminate accumulate
  bool done;
  do {
    done = sender->test();
    done &= accumer->test();
#ifndef USE_SERVER_THREAD
    // in case no thread is running behind, we need to cycle this to flush
    size_t d = done ? 0 : 1;
    mpi__->soft_allreduce(&d, 1);
    done = d == 0;
    if (!done) sender->flush();
    if (!done) accumer->flush();
#endif
    if (!done) this_thread::sleep_for(sleeptime__);
  } while (!done);
}


void DistMatrix::rotate(const int i, const int j, const double gamma) {
  rotate(vector<tuple<int, int, double>>(1, make_tuple(i, j, gamma)));
}


double DistMatrix::ddot(const DistMatrix& o) const {
  assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
  double sum = size() ? ddot_(size(), local_.get(), 1, o.local_.get(), 1) : 0.0;
  mpi__->allreduce(&sum, 1);
  return sum;
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

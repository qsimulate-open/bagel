//
// BAGEL - Parallel electron correlation program.
// Filename: distcivec.cc
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

// TODO until GCC fixes this bug
#define _GLIBCXX_USE_NANOSLEEP

#include <sstream>
#include <src/fci/civec.h>
#include <src/math/algo.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

DistCivec::DistCivec(shared_ptr<const Determinants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()), dist_(lena_, mpi__->size()) {
  tie(astart_, aend_) = dist_.range(mpi__->rank());

  alloc_ = size();
  local_ = unique_ptr<double[]>(new double[alloc_]);
  fill_n(local_.get(), alloc_, 0.0);

  mutex_ = vector<mutex>(asize());
}


DistCivec::DistCivec(const DistCivec& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_), dist_(lena_, mpi__->size()) {
  tie(astart_, aend_) = dist_.range(mpi__->rank());
  alloc_ = size();
  local_ = unique_ptr<double[]>(new double[alloc_]);
  copy_n(o.local_.get(), alloc_, local_.get());

  mutex_ = vector<mutex>(asize());
}


DistCivec& DistCivec::operator=(const DistCivec& o) {
  assert(o.size() == size());
  copy_n(o.local_.get(), alloc_, local_.get());
  return *this;
}


void DistCivec::init_mpi_accumulate() const {
  send_  = make_shared<SendRequest>();
  accum_ = make_shared<AccRequest>(local_.get(), &mutex_);
}


void DistCivec::init_mpi_recv() const {
  put_   = make_shared<PutRequest>(local_.get());
  recv_  = make_shared<RecvRequest>();
}


void DistCivec::accumulate_bstring_buf(unique_ptr<double[]>& buf, const size_t a) const {
  assert(accum_ && send_);
  const size_t mpirank = mpi__->rank();
  size_t rank, off;
  tie(rank, off) = dist_.locate(a);

  if (mpirank == rank) {
    lock_guard<mutex> lock(mutex_[off]);
    daxpy_(lenb_, 1.0, buf.get(), 1, local_.get()+off*lenb_, 1);
  } else {
    send_->request_send(move(buf), lenb_, rank, off*lenb_);
  }
}


int DistCivec::get_bstring_buf(double* buf, const size_t a) const {
  assert(put_ && recv_);
  const size_t mpirank = mpi__->rank();
  size_t rank, off;
  tie(rank, off) = dist_.locate(a);

  int out = -1;
  if (mpirank == rank) {
    copy_n(local_.get()+off*lenb_, lenb_, buf);
  } else {
    out = recv_->request_recv(buf, lenb_, rank, off*lenb_);
  }
  return out;
}


void DistCivec::terminate_mpi_accumulate() const {
  assert(accum_ && send_);

  bool done;
  do {
    done = send_->test();
    done &= accum_->test();
#ifndef USE_SERVER_THREAD
    // in case no thread is running behind, we need to cycle this to flush
    size_t d = done ? 0 : 1;
    mpi__->soft_allreduce(&d, 1);
    done = d == 0;
    if (!done) send_->flush();
    if (!done) accum_->flush();
#endif
    if (!done) this_thread::sleep_for(sleeptime__);
  } while (!done);

  // cancel all MPI calls
  send_  = shared_ptr<SendRequest>();
  accum_ = shared_ptr<AccRequest>();
}


void DistCivec::terminate_mpi_recv() const {
  assert(put_ && recv_);

  bool done;
  do {
    done = recv_->test();
#ifndef USE_SERVER_THREAD
    // in case no thread is running behind, we need to cycle this to flush
    size_t d = done ? 0 : 1;
    mpi__->soft_allreduce(&d, 1);
    done = d == 0;
    if (!done) put_->flush();
#endif
    if (!done) this_thread::sleep_for(sleeptime__);
  } while (!done);

  // cancel all MPI calls
  recv_  = shared_ptr<RecvRequest>();
  put_   = shared_ptr<PutRequest>();
}


double DistCivec::dot_product(const DistCivec& o) const {
  assert(size() == o.size());
  double sum = size() ? ddot_(size(), local(), 1, o.local(), 1) : 0.0;
  mpi__->allreduce(&sum, 1);
  return sum;
}


double DistCivec::norm() const {
  return sqrt(dot_product(*this));
}


void DistCivec::daxpy(const double a, const DistCivec& o) {
  assert(size() == o.size());
  for (size_t i = 0; i != asize(); ++i) {
    lock_guard<mutex> lock(mutex_[i]);
    daxpy_(lenb_, a, o.local()+i*lenb_, 1, local()+i*lenb_, 1);
  }
}


void DistCivec::scale(const double a) {
  for (size_t i = 0; i != asize(); ++i) {
    lock_guard<mutex> lock(mutex_[i]);
    dscal_(lenb_, a, local()+i*lenb_, 1);
  }
}


double DistCivec::orthog(list<shared_ptr<const DistCivec>> c) {
  for (auto& iter : c)
    project_out(iter);
  const double norm = this->norm();
  const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
  scale(scal);
  return norm;
}


double DistCivec::orthog(shared_ptr<const DistCivec> o) {
  list<shared_ptr<const DistCivec>> v = {o};
  return orthog(v);
}


double DistCivec::variance() const {
  return dot_product(*this) / (lena_*lenb_);
}


shared_ptr<DistCivec> DistCivec::transpose() const {
  shared_ptr<Determinants> det = det_->transpose();
  auto out = make_shared<DistCivec>(det);

  const size_t myrank = mpi__->rank();

  // transpose each segment
  shared_ptr<DistCivec> trans = clone();
  for (int i = 0; i != mpi__->size(); ++i) {
    tuple<size_t, size_t> outrange = out->dist_.range(i);
    tuple<size_t, size_t> thisrange = dist_.range(i);

    unique_ptr<double[]> tmp(new double[out->dist_.size(i)*asize()]);
    for (size_t j = 0; j != asize(); ++j)
      copy_n(local()+get<0>(outrange)+j*lenb_, out->dist_.size(i), tmp.get()+j*out->dist_.size(i));

    const size_t off = get<0>(outrange)*asize();
    copy_n(tmp.get(), out->dist_.size(i)*asize(), trans->local()+off);
    if (i != myrank) {
      out->transp_.push_back(mpi__->request_send(trans->local()+off, out->dist_.size(i)*asize(), i, myrank));
      out->transp_.push_back(mpi__->request_recv(out->local()+out->asize()*get<0>(thisrange), out->asize()*dist_.size(i), i, i));
    } else {
      copy_n(trans->local()+off, out->asize()*asize(), out->local()+astart()*out->asize());
    }
  }

  // keep trans
  out->buf_ = trans;
  return out;
}


void DistCivec::transpose_wait() {
  for (auto& i : transp_)
    mpi__->wait(i);

  buf_ = shared_ptr<DistCivec>();
  buf_ = clone();

  assert(buf_);
  mytranspose_(local(), asize(), lenb_, buf_->local());
  copy_n(buf_->local(), asize()*lenb_, local());

  buf_ = shared_ptr<DistCivec>();
}

shared_ptr<Civector<double>> DistCivec::civec() const { return make_shared<Civector<double>>(*this); }

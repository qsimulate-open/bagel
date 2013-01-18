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


#include <sstream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp> 
#include <src/fci/civec.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

DistCivec::DistCivec(shared_ptr<const Determinants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()), win_(-1), dist_(lena_, mpi__->size()) {
  tie(astart_, aend_) = dist_.range(mpi__->rank());

  alloc_ = size();
  local_ = unique_ptr<double[]>(new double[alloc_]);
  fill_n(local_.get(), alloc_, 0.0);

  mutex_ = vector<boost::mutex>(asize());
}


DistCivec::DistCivec(const DistCivec& o) : det_(o.det_), lena_(o.lena_), lenb_(o.lenb_), win_(-1), dist_(lena_, mpi__->size()) {
  tie(astart_, aend_) = dist_.range(mpi__->rank());
  alloc_ = size();
  local_ = unique_ptr<double[]>(new double[alloc_]);
  copy_n(o.local_.get(), alloc_, local_.get());

  mutex_ = vector<boost::mutex>(asize());
}


DistCivec& DistCivec::operator=(const DistCivec& o) {
  assert(o.size() == size());
  copy_n(o.local_.get(), alloc_, local_.get());
  return *this;
} 


#if 0
void DistCivec::open_window() const {
  assert(win_ == -1);
  win_ = mpi__->win_create(local_.get(), alloc_);
}


void DistCivec::close_window() const {
  assert(win_ != -1);
  mpi__->win_free(win_);
  win_ = -1; 
}


void DistCivec::fence() const {
  assert(win_ != -1);
  mpi__->win_fence(win_);
}


void DistCivec::get_bstring(double* buf, const size_t a) const {
  const size_t mpirank = mpi__->rank();
  size_t rank, off;
  tie(rank, off) = dist_.locate(a);

  assert(win_ != -1);
  if (mpirank == rank) {
    copy_n(local_.get()+off*lenb_, lenb_, buf);
  } else {
    mpi__->get(buf, lenb_, rank, off*lenb_, win_); 
  }
}
#endif


void DistCivec::init_mpi_accumulate() const {
  send_  = shared_ptr<SendRequest>(new SendRequest());
  accum_ = shared_ptr<AccRequest>(new AccRequest(local_.get(), &mutex_));
  accum_->init_request();
}


void DistCivec::init_mpi_recv() const {
  put_   = shared_ptr<PutRequest>(new PutRequest(local_.get()));
  recv_  = shared_ptr<RecvRequest>(new RecvRequest());
  put_->init_request();
}


void DistCivec::accumulate_bstring_buf(unique_ptr<double[]>& buf, const size_t a) const {
  assert(accum_ && send_);
  const size_t mpirank = mpi__->rank();
  size_t rank, off;
  tie(rank, off) = dist_.locate(a);

  if (mpirank == rank) {
    boost::lock_guard<boost::mutex> lock(mutex_[off]);
    daxpy_(lenb_, 1.0, buf.get(), 1, local_.get()+off*lenb_, 1);
  } else {
    send_->request_send(std::move(buf), lenb_, rank, off*lenb_);
  }
}


void DistCivec::get_bstring_buf(double* buf, const size_t a) const {
  assert(put_ && recv_);
  const size_t mpirank = mpi__->rank();
  size_t rank, off;
  tie(rank, off) = dist_.locate(a);

  if (mpirank == rank) {
    copy_n(local_.get()+off*lenb_, lenb_, buf);
  } else {
    recv_->request_recv(buf, lenb_, rank, off*lenb_);
  }
}


void DistCivec::flush_accumulate() const {
  assert(accum_ && send_);
  send_->flush();
  accum_->flush(); 
}


void DistCivec::flush_recv() const {
  assert(put_);
//put_->flush();
}


void DistCivec::recv_wait() const {
  assert(put_ && recv_);
  bool done;
  do {
    done = recv_->test1();
    done &= recv_->test2();
    if (!done) boost::this_thread::sleep(boost::posix_time::milliseconds(1));
  } while (!done);
}


void DistCivec::terminate_mpi_accumulate() const {
  assert(accum_ && send_);

  bool done;
  do {
    done = send_->test1();
    accum_->flush();
    done &= send_->test2();
    done &= accum_->test2();
    done &= send_->test3();
    done &= accum_->test3();
    int d = done ? 0 : 1; 
    mpi__->allreduce(&d, 1);
    done = d == 0;
    if (!done) boost::this_thread::sleep(boost::posix_time::milliseconds(1));
  } while (!done);

  // cancel all MPI calls
  send_  = shared_ptr<SendRequest>();
  accum_ = shared_ptr<AccRequest>();
}


void DistCivec::terminate_mpi_recv() const {
  assert(put_ && recv_);

  bool done;
  do {
    done = recv_->test1();
    done &= recv_->test2();
    if (!done) boost::this_thread::sleep(boost::posix_time::milliseconds(1));
  } while (!done);

  // cancel all MPI calls
  recv_  = shared_ptr<RecvRequest>();
  put_   = shared_ptr<PutRequest>();
}


double DistCivec::ddot(const DistCivec& o) const {
  assert(size() == o.size());
  double sum = size() ? ddot_(size(), local_.get(), 1, o.local_.get(), 1) : 0.0;
  mpi__->allreduce(&sum, 1);
  return sum;
}


double DistCivec::norm() const {
  return std::sqrt(ddot(*this));
}


void DistCivec::daxpy(const double a, const DistCivec& o) {
  assert(size() == o.size());
  daxpy_(size(), a, o.local(), 1, local(), 1);
}


void DistCivec::scale(const double a) {
  dscal_(size(), a, local(), 1);
}


double DistCivec::orthog(list<shared_ptr<const DistCivec> > c) {
  for (auto& iter : c)
    project_out(iter);
  const double norm = this->norm();
  const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
  scale(scal);
  return 1.0/scal;
}


double DistCivec::orthog(shared_ptr<const DistCivec> o) {
  list<shared_ptr<const DistCivec> > v = {o};
  return orthog(v);
}


double DistCivec::variance() const {
  return ddot(*this) / (lena_*lenb_);
}

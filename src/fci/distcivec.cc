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
#include <src/fci/civec.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

DistCivec::DistCivec(shared_ptr<const Determinants> det) : det_(det), lena_(det->lena()), lenb_(det->lenb()), win_(-1), dist_(lena_, mpi__->size()) {
  tie(astart_, aend_) = dist_.range(mpi__->rank()); 

  alloc_ = size()*lenb_;
  local_ = unique_ptr<double[]>(new double[alloc_]);
  fill_n(local_.get(), alloc_, 0.0);

  accum_ = shared_ptr<AccRequest>(new AccRequest(0, local_.get()));
  send_  = shared_ptr<SendRequest>(new SendRequest());
}


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


void DistCivec::put_bstring(const double* buf, const size_t a) const {
  const size_t mpirank = mpi__->rank();
  size_t rank, off;
  tie(rank, off) = dist_.locate(a);

  assert(win_ != -1);
  if (mpirank == rank) {
    copy_n(buf, lenb_, local_.get()+off*lenb_);
  } else {
    mpi__->put(buf, lenb_, rank, off*lenb_, win_); 
  }
}


void DistCivec::accumulate_bstring(const double* buf, const size_t a) const {
  const size_t mpirank = mpi__->rank();
  size_t rank, off;
  tie(rank, off) = dist_.locate(a);

  assert(win_ != -1);
  if (mpirank == rank) {
    daxpy_(lenb_, 1.0, buf, 1, local_.get()+off*lenb_, 1);
  } else {
    mpi__->accumulate(buf, lenb_, rank, off*lenb_, win_); 
  }
}


void DistCivec::accumulate_bstring_buf(unique_ptr<double[]>& buf, const size_t a) const {
  const size_t mpirank = mpi__->rank();
  size_t rank, off;
  tie(rank, off) = dist_.locate(a);

  if (mpirank == rank) {
    daxpy_(lenb_, 1.0, buf.get(), 1, local_.get()+off*lenb_, 1);
  } else {
    send_->request_send(std::move(buf), lenb_, rank, off); 
  }
}


void DistCivec::flush() const {
  send_->flush();
  accum_->flush(); 
}


void DistCivec::wait() const {
  send_->wait();
  accum_->wait();
}

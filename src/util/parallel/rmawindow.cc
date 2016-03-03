//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mrawindow.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <cassert>
#include <stdexcept>
#include <src/util/math/algo.h>
#include <src/util/parallel/rmawindow.h>
#include <src/util/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;


template<typename DataType>
void RMATask<DataType>::wait() {
#ifdef HAVE_MPI_H
  MPI_Wait(&tag, MPI_STATUS_IGNORE);
#endif
}


template<typename DataType>
bool RMATask<DataType>::test() {
  int flag = 0;
#ifdef HAVE_MPI_H
  MPI_Test(&tag, &flag, MPI_STATUS_IGNORE);
#endif
  return flag;
}


template<typename DataType>
RMAWindow<DataType>::RMAWindow() : initialized_(false) {
#ifndef HAVE_MPI_H
  throw logic_error("RMAWindow should be used with MPI");
#endif
}


template<typename DataType>
void RMAWindow<DataType>::initialize() {
#ifdef HAVE_MPI_H
  assert(!initialized_);
  // allocate a window
  MPI_Win_allocate(localsize()*sizeof(DataType), sizeof(DataType), MPI_INFO_NULL, MPI_COMM_WORLD, &win_base_, &win_);
  MPI_Win_lock_all(MPI_MODE_NOCHECK, win_);

  initialized_ = true;
  zero();
#endif
}


template<typename DataType>
RMAWindow<DataType>::~RMAWindow() {
#ifdef HAVE_MPI_H
  assert(initialized_);
  MPI_Win_unlock_all(win_);
  MPI_Win_free(&win_);
#endif
}


template<typename DataType>
void RMAWindow<DataType>::zero() {
  assert(initialized_);
  fence();
  const size_t loc = localsize();
  if (loc)
    fill_n(win_base_, loc, 0.0);
  fence_local();
  mpi__->barrier();
}


template<typename DataType>
void RMAWindow<DataType>::scale(const DataType& a) {
  assert(initialized_);
  fence();
  const size_t loc = localsize();
  if (loc)
    blas::scale_n(a, win_base_, loc);
  fence_local();
  mpi__->barrier();
}


template<typename DataType>
RMAWindow<DataType>& RMAWindow<DataType>::operator=(const RMAWindow<DataType>& o) {
  assert(o.initialized_);
  if (!initialized_ && o.initialized_)
    initialize();
  fence();
  o.fence();
  const size_t loc = localsize();
  assert(loc == o.localsize());
  if (loc)
    copy_n(o.win_base_, loc, win_base_);
  fence_local();
  o.fence_local();
  mpi__->barrier();
  return *this;
}


template<typename DataType>
void RMAWindow<DataType>::fence() const {
#ifdef HAVE_MPI_H
  assert(initialized_);
  MPI_Win_flush_all(win_);
  mpi__->barrier();
#endif
}


template<typename DataType>
void RMAWindow<DataType>::fence_local() const {
#ifdef HAVE_MPI_H
  assert(initialized_);
  MPI_Win_sync(win_);
#endif
}


template<typename DataType>
void RMAWindow<DataType>::ax_plus_y(const DataType& a, const RMAWindow<DataType>& o) {
  assert(initialized_);
  fence();
  o.fence();
  const size_t loc = localsize();
  if (loc)
    blas::ax_plus_y_n(a, o.win_base_, loc, win_base_);
  fence_local();
  o.fence_local();
  mpi__->barrier();
}


template<typename DataType>
DataType RMAWindow<DataType>::dot_product(const RMAWindow<DataType>& o) const {
  assert(initialized_);
  fence();
  o.fence();
  const size_t loc = localsize();
  DataType out = loc ? blas::dot_product(win_base_, loc, o.win_base_) : 0.0;
  fence_local();
  o.fence_local();
  mpi__->allreduce(&out, 1);
  return out;
}


template<typename DataType>
unique_ptr<DataType[]> RMAWindow<DataType>::rma_get(const size_t& key) const {
  assert(initialized_);
#ifdef HAVE_MPI_H
  size_t rank, off, size;
  tie(rank, off, size) = locate(key);
  unique_ptr<DataType[]> out(new DataType[size]);
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  MPI_Get(out.get(), size, type, rank, off, size, type, win_);
  MPI_Win_flush_local(rank, win_);
  return move(out);
#else
  return nullptr;
#endif
}


template<typename DataType>
void RMAWindow<DataType>::rma_put(const unique_ptr<DataType[]>& dat, const size_t& key) {
  assert(initialized_);
#ifdef HAVE_MPI_H
  size_t rank, off, size;
  tie(rank, off, size) = locate(key);
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  MPI_Put(dat.get(), size, type, rank, off, size, type, win_);
  MPI_Win_flush_local(rank, win_);
#endif
}


template<typename DataType>
void RMAWindow<DataType>::rma_add(const unique_ptr<DataType[]>& dat, const size_t& key) {
  assert(initialized_);
#ifdef HAVE_MPI_H
  size_t rank, off, size;
  tie(rank, off, size) = locate(key);
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  MPI_Accumulate(dat.get(), size, type, rank, off, size, type, MPI_SUM, win_);
  MPI_Win_flush_local(rank, win_);
#endif
}


template<typename DataType>
shared_ptr<RMATask<DataType>> RMAWindow<DataType>::rma_rget(DataType* buf, const size_t& key) const {
  shared_ptr<RMATask<DataType>> out;
#ifdef HAVE_MPI_H
  size_t rank, off, size;
  tie(rank, off, size) = locate(key);
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  MPI_Request req;
  MPI_Rget(buf, size, type, rank, off, size, type, win_, &req);
  out = make_shared<RMATask<DataType>>(move(req));
#endif
  return out;

}


template<typename DataType>
shared_ptr<RMATask<DataType>> RMAWindow<DataType>::rma_radd(unique_ptr<DataType[]>&& buf, const size_t& key) {
  shared_ptr<RMATask<DataType>> out;
#ifdef HAVE_MPI_H
  size_t rank, off, size;
  tie(rank, off, size) = locate(key);
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  MPI_Request req;
  MPI_Raccumulate(buf.get(), size, type, rank, off, size, type, MPI_SUM, win_, &req);
  out = make_shared<RMATask<DataType>>(move(req), move(buf));
#endif
  return out;
}


template<typename DataType>
void RMAWindow<DataType>::set_element(const size_t rank, const size_t disp, const DataType a) {
#ifdef HAVE_MPI_H
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  MPI_Put(&a, 1, type, rank, disp, 1, type, win_);
  MPI_Win_flush_local(rank, win_);
#endif
}


template<typename DataType>
void RMAWindow<DataType>::accumulate_buffer(const DataType a, const unique_ptr<DataType[]>& buf) {
  fence();
  blas::ax_plus_y_n(a, buf.get(), localsize(), win_base_);
  fence_local();
  mpi__->barrier();
}


template class RMATask<double>;
template class RMATask<complex<double>>;
template class RMAWindow<double>;
template class RMAWindow<complex<double>>;

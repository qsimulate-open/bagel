//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: storage.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <ga.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <src/util/f77.h>
#include <src/util/math/algo.h>
#include <src/smith/storage.h>
#include <src/util/parallel/mpi_interface.h>

using namespace bagel::SMITH;
using namespace std;

// size contains hashkey and length (in this order)
template<typename DataType>
StorageIncore<DataType>::StorageIncore(const map<size_t, size_t>& size, bool init) : initialized_(false) {
  static_assert(is_same<DataType, double>::value or is_same<DataType, complex<double>>::value, "illegal Type in StorageIncore");

  // first prepare some variables
  totalsize_ = 0;
  for (auto& i : size) {
    if (i.second > 0)
      hashtable_.emplace(i.first, make_pair(totalsize_, totalsize_+i.second));
    totalsize_ += i.second;
  }

  // store block data
  const size_t blocksize = (totalsize_-1)/mpi__->size()+1;
  size_t tsize = 0;
  for (auto& i : size) {
    if (i.second == 0) continue;
    if (blocks_.size()*blocksize <= tsize)
      blocks_.emplace(tsize, blocks_.size());
    tsize += i.second;
  }
  assert(totalsize_ == tsize);

  // set local_lo_ and hi_
  local_lo_ = local_hi_ = numeric_limits<size_t>::max();
  for (auto iter = blocks_.begin(); iter != blocks_.end(); ++iter)
    if (iter->second == mpi__->rank()) {
      local_lo_ = iter->first;
      auto iterp = ++iter;
      local_hi_ = iterp != blocks_.end() ? iterp->first : totalsize_;
      break;
    }

  // here we initialize the global array storage
  if (init)
    initialize();
}


template<typename DataType>
void StorageIncore<DataType>::initialize() {
  assert(!initialized_);
  // allocate a window
  MPI_Aint size = localsize()*sizeof(DataType);
  MPI_Win_allocate(size, sizeof(DataType), MPI_INFO_NULL, MPI_COMM_WORLD, &win_base_, &win_);

  initialized_ = true;
  zero();
}


template<typename DataType>
StorageIncore<DataType>::~StorageIncore() {
  MPI_Win_fence(0, win_);
  MPI_Win_free(&win_);
}


template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block_(const size_t& key) const {
  assert(initialized_);
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  size_t rank, off, size;
  tie(rank, off, size) = locate(key);

  unique_ptr<DataType[]> out(new DataType[size]);
  if (rank != mpi__->rank())
    MPI_Get(out.get(), size, type, rank, off, size, type, win_);
  else
    copy_n(win_base_+off, size, out.get());
  return move(out);
}


template<typename DataType>
void StorageIncore<DataType>::put_block_(const unique_ptr<DataType[]>& dat, const size_t& key) {
  assert(initialized_);
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  size_t rank, off, size;
  tie(rank, off, size) = locate(key);
  if (rank != mpi__->rank())
    MPI_Put(dat.get(), size, type, rank, off, size, type, win_);
  else
    copy_n(dat.get(), size, win_base_+off);
}


template<typename DataType>
void StorageIncore<DataType>::add_block_(const unique_ptr<DataType[]>& dat, const size_t& key) {
  assert(initialized_);
  auto type = is_same<double,DataType>::value ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX;
  size_t rank, off, size;
  tie(rank, off, size) = locate(key);

  if (rank != mpi__->rank())
    MPI_Accumulate(dat.get(), size, type, rank, off, size, type, MPI_SUM, win_);
  else
    blas::ax_plus_y_n(1.0, dat.get(), size, win_base_+off);
}


template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block() const {
  return get_block_(generate_hash_key());
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0) const {
  return get_block_(generate_hash_key(i0));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1) const {
  return get_block_(generate_hash_key(i0, i1));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2) const {
  return get_block_(generate_hash_key(i0, i1, i2));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3) const {
  return get_block_(generate_hash_key(i0, i1, i2, i3));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                          const Index& i4) const {
  return get_block_(generate_hash_key(i0, i1, i2, i3, i4));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                          const Index& i4, const Index& i5) const {
  return get_block_(generate_hash_key(i0, i1, i2, i3, i4, i5));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                          const Index& i4, const Index& i5, const Index& i6) const {
  return get_block_(generate_hash_key(i0, i1, i2, i3, i4, i5, i6));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                          const Index& i4, const Index& i5, const Index& i6, const Index& i7) const {
  return get_block_(generate_hash_key(i0, i1, i2, i3, i4, i5, i6, i7));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(vector<Index> i) const {
  return get_block_(generate_hash_key(i));
}


template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat) {
  put_block_(dat, generate_hash_key());
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0) {
  put_block_(dat, generate_hash_key(i0));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) {
  put_block_(dat, generate_hash_key(i0, i1));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) {
  put_block_(dat, generate_hash_key(i0, i1, i2));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3, i4));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5, const Index& i6) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6, i7));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, vector<Index> i) {
  put_block_(dat, generate_hash_key(i));
}


template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat) {
  add_block_(dat, generate_hash_key());
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0) {
  add_block_(dat, generate_hash_key(i0));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) {
  add_block_(dat, generate_hash_key(i0, i1));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) {
  add_block_(dat, generate_hash_key(i0, i1, i2));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  add_block_(dat, generate_hash_key(i0, i1, i2, i3));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4) {
  add_block_(dat, generate_hash_key(i0, i1, i2, i3, i4));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5) {
  add_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5, const Index& i6) {
  add_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  add_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6, i7));
}


template<typename DataType>
void StorageIncore<DataType>::zero() {
  assert(initialized_);
  const size_t loc = localsize();
  if (loc)
    fill_n(win_base_, loc, 0.0);
  MPI_Win_fence(0, win_);
}


template<typename DataType>
void StorageIncore<DataType>::scale(const DataType& a) {
  assert(initialized_);
  const size_t loc = localsize();
  if (loc)
    blas::scale_n(a, win_base_, loc);
  MPI_Win_fence(0, win_);
}


template<typename DataType>
StorageIncore<DataType>& StorageIncore<DataType>::operator=(const StorageIncore<DataType>& o) {
  if (!initialized_ && o.initialized_)
    initialize();
  const size_t loc = localsize();
  if (loc)
    copy_n(o.win_base_, loc, win_base_);
  MPI_Win_fence(0, win_);
  return *this;
}


template<typename DataType>
tuple<size_t,size_t,size_t> StorageIncore<DataType>::locate(const size_t key) const {
  auto iter = hashtable_.find(key);
  if (iter == hashtable_.end())
    throw logic_error("a key was not found in Storage::locate(const size_t key)");
  const size_t lo = iter->second.first;
  const size_t hi = iter->second.second;
  auto b = blocks_.upper_bound(lo);
  --b;
  return make_tuple(b->second, lo - b->first, hi - lo);
}


template<typename DataType>
bool StorageIncore<DataType>::is_local(const size_t key) const {
  auto iter = hashtable_.find(key);
  assert(iter != hashtable_.end());
  const size_t lo = iter->second.first;
  return lo >= local_lo_ && lo < local_hi_;
}


template<typename DataType>
size_t StorageIncore<DataType>::localsize() const {
  return local_hi_ - local_lo_;
}


template<typename DataType>
void StorageIncore<DataType>::ax_plus_y(const DataType& a, const StorageIncore<DataType>& o) {
  assert(initialized_);
  const size_t loc = localsize();
  if (loc)
    blas::ax_plus_y_n(a, o.win_base_, loc, win_base_);
  MPI_Win_fence(0, win_);
}


template<typename DataType>
DataType StorageIncore<DataType>::dot_product(const StorageIncore<DataType>& o) const {
  assert(initialized_);
  MPI_Win_fence(0, win_);
  const size_t loc = localsize();
  DataType out = loc ? blas::dot_product(win_base_, loc, o.win_base_) : 0.0;
  mpi__->allreduce(&out, 1);
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class StorageIncore<double>;
template class StorageIncore<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

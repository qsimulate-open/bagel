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
StorageIncore<DataType>::StorageIncore(const map<size_t, size_t>& size, bool init) : RMAWindow<DataType>() {
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
unique_ptr<DataType[]> StorageIncore<DataType>::get_block() const {
  return rma_get(generate_hash_key());
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0) const {
  return rma_get(generate_hash_key(i0));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1) const {
  return rma_get(generate_hash_key(i0, i1));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2) const {
  return rma_get(generate_hash_key(i0, i1, i2));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3) const {
  return rma_get(generate_hash_key(i0, i1, i2, i3));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                          const Index& i4) const {
  return rma_get(generate_hash_key(i0, i1, i2, i3, i4));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                          const Index& i4, const Index& i5) const {
  return rma_get(generate_hash_key(i0, i1, i2, i3, i4, i5));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                          const Index& i4, const Index& i5, const Index& i6) const {
  return rma_get(generate_hash_key(i0, i1, i2, i3, i4, i5, i6));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                          const Index& i4, const Index& i5, const Index& i6, const Index& i7) const {
  return rma_get(generate_hash_key(i0, i1, i2, i3, i4, i5, i6, i7));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(vector<Index> i) const {
  return rma_get(generate_hash_key(i));
}


template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat) {
  rma_put(dat, generate_hash_key());
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0) {
  rma_put(dat, generate_hash_key(i0));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) {
  rma_put(dat, generate_hash_key(i0, i1));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) {
  rma_put(dat, generate_hash_key(i0, i1, i2));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  rma_put(dat, generate_hash_key(i0, i1, i2, i3));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4) {
  rma_put(dat, generate_hash_key(i0, i1, i2, i3, i4));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5) {
  rma_put(dat, generate_hash_key(i0, i1, i2, i3, i4, i5));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5, const Index& i6) {
  rma_put(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  rma_put(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6, i7));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(const unique_ptr<DataType[]>& dat, vector<Index> i) {
  rma_put(dat, generate_hash_key(i));
}


template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat) {
  rma_add(dat, generate_hash_key());
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0) {
  rma_add(dat, generate_hash_key(i0));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) {
  rma_add(dat, generate_hash_key(i0, i1));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) {
  rma_add(dat, generate_hash_key(i0, i1, i2));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  rma_add(dat, generate_hash_key(i0, i1, i2, i3));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4) {
  rma_add(dat, generate_hash_key(i0, i1, i2, i3, i4));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5) {
  rma_add(dat, generate_hash_key(i0, i1, i2, i3, i4, i5));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5, const Index& i6) {
  rma_add(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6));
}

template<typename DataType>
void StorageIncore<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  rma_add(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6, i7));
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class bagel::SMITH::StorageIncore<double>;
template class bagel::SMITH::StorageIncore<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_CLASS_EXPORT_IMPLEMENT(StorageIncore<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(StorageIncore<complex<double>>)

#endif

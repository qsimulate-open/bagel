//
// BAGEL - Parallel electron correlation program.
// Filename: storage.cc
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/util/f77.h>
#include <src/util/math/algo.h>
#include <src/smith/storage.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace bagel::SMITH;
using namespace std;


template<typename DataType>
StorageBlock<DataType>::StorageBlock(const size_t size, const bool init) : size_(size), initialized_(init) {
  static_assert(is_same<DataType, double>::value or is_same<DataType, complex<double>>::value, "illegal Type in StorageBlock");
  if (init) {
    data_ = unique_ptr<DataType[]>(new DataType[size_]);
    zero();
  }
}


template<typename DataType>
void StorageBlock<DataType>::zero() {
  if (initialized_)
    fill_n(data(), size_, 0.0);
}


template<typename DataType>
StorageBlock<DataType>& StorageBlock<DataType>::operator=(const StorageBlock<DataType>& o) {
  if (o.initialized_ && !initialized_) {
    data_ = unique_ptr<DataType[]>(new DataType[size_]);
    initialized_ = true;
  }
  if (o.initialized_)
    copy_n(o.data(), size_, data());
  return *this;
}


template<typename DataType>
void StorageBlock<DataType>::put_block(unique_ptr<DataType[]>&& o) {
  assert(!initialized_);
  initialized_ = true;
  data_ = move(o);
}


template<typename DataType>
void StorageBlock<DataType>::add_block(const std::unique_ptr<DataType[]>& o) {
  assert(initialized_);
  blas::ax_plus_y_n(1.0, o.get(), size_, data());
}


template<typename DataType>
unique_ptr<DataType[]> StorageBlock<DataType>::get_block() const {
  assert(initialized_);
  unique_ptr<DataType[]> out(new DataType[size_]);
  copy_n(data_.get(), size_, out.get());
  return move(out);
}


template<typename DataType>
unique_ptr<DataType[]> StorageBlock<DataType>::move_block() {
  if (!initialized_) {
    initialized_ = true;
    data_ = unique_ptr<DataType[]>(new DataType[size_]);
    zero();
  }
  initialized_ = false;
  return move(data_);
}


template<typename DataType>
DataType StorageBlock<DataType>::dot_product(const StorageBlock<DataType>& o) const {
  assert(size_ == o.size_ && !(initialized_ ^ o.initialized_));
  return initialized_ ? blas::dot_product(data(), size_, o.data()) : DataType(0.0);
}


template<typename DataType>
void StorageBlock<DataType>::ax_plus_y(const DataType& a, const StorageBlock<DataType>& o) {
  assert(size_ == o.size_ && !(initialized_ ^ o.initialized_));
  if (initialized_)
    blas::ax_plus_y_n(a, o.data(), size_, data());
}


template<typename DataType>
void StorageBlock<DataType>::scale(const DataType& a) {
  if (initialized_)
    blas::scale_n(a, data(), size_);
}


template<typename DataType>
void StorageBlock<DataType>::conjugate_inplace() {
  blas::conj_n(data_.get(), size_alloc());
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename DataType>
StorageIncore<DataType>::StorageIncore(const map<size_t, size_t>& size, bool init) : Storage_base<StorageBlock<DataType>>(size, init) {
  static_assert(is_same<DataType, double>::value or is_same<DataType, complex<double>>::value, "illegal Type in StorageIncore");
}


template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::get_block_(const size_t& key) const {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::get_block(const size_t&)");
  return hash->second->get_block();
}


template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block_(const size_t& key) {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::move_block(const size_t&)");
  return hash->second->move_block();
}


template<typename DataType>
void StorageIncore<DataType>::put_block_(unique_ptr<DataType[]>& dat, const size_t& key) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");
  hash->second->put_block(move(dat));
}


template<typename DataType>
void StorageIncore<DataType>::add_block_(const unique_ptr<DataType[]>& dat, const size_t& key) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");
  hash->second->add_block(dat);
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
unique_ptr<DataType[]> StorageIncore<DataType>::get_block(vector<Index> i) {
  return get_block_(generate_hash_key(i));
}


template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block() {
  return move_block_(generate_hash_key());
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block(const Index& i0) {
  return move_block_(generate_hash_key(i0));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block(const Index& i0, const Index& i1) {
  return move_block_(generate_hash_key(i0, i1));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2) {
  return move_block_(generate_hash_key(i0, i1, i2));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  return move_block_(generate_hash_key(i0, i1, i2, i3));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4) {
  return move_block_(generate_hash_key(i0, i1, i2, i3, i4));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5) {
  return move_block_(generate_hash_key(i0, i1, i2, i3, i4, i5));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5, const Index& i6) {
  return move_block_(generate_hash_key(i0, i1, i2, i3, i4, i5, i6));
}

template<typename DataType>
unique_ptr<DataType[]> StorageIncore<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  return move_block_(generate_hash_key(i0, i1, i2, i3, i4, i5, i6, i7));
}


template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat) {
  put_block_(dat, generate_hash_key());
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0) {
  put_block_(dat, generate_hash_key(i0));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) {
  put_block_(dat, generate_hash_key(i0, i1));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) {
  put_block_(dat, generate_hash_key(i0, i1, i2));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                     const Index& i4) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3, i4));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                     const Index& i4, const Index& i5) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                     const Index& i4, const Index& i5, const Index& i6) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                     const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  put_block_(dat, generate_hash_key(i0, i1, i2, i3, i4, i5, i6, i7));
}

template<typename DataType>
void StorageIncore<DataType>::put_block(unique_ptr<DataType[]>& dat, vector<Index> i) {
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
  for (auto& i : hashtable_)
    i.second->zero();
}


template<typename DataType>
void StorageIncore<DataType>::scale(const DataType& a) {
  for (auto& i : hashtable_)
    i.second->scale(a);
}


template<typename DataType>
StorageIncore<DataType>& StorageIncore<DataType>::operator=(const StorageIncore<DataType>& o) {
  if (hashtable_.size() == o.hashtable_.size()) {
    auto i = hashtable_.begin();
    for (auto& j : o.hashtable_) {
      *i->second = *j.second;
      ++i;
    }
  } else {
    throw logic_error("Trying to copy something different in StorageIncore");
  }
  return *this;
}


template<typename DataType>
void StorageIncore<DataType>::ax_plus_y(const DataType& a, const StorageIncore<DataType>& o) {
  if (hashtable_.size() == o.hashtable_.size()) {
    auto i = hashtable_.begin();
    for (auto& j : o.hashtable_) {
      i->second->ax_plus_y(a, *j.second);
      ++i;
    }
  } else {
    throw logic_error("Trying to copy something different in StorageIncore");
  }
}


template<typename DataType>
DataType StorageIncore<DataType>::dot_product(const StorageIncore<DataType>& o) const {
  DataType out = 0.0;
  if (hashtable_.size() == o.hashtable_.size()) {
    auto i = hashtable_.begin();
    for (auto& j : o.hashtable_) {
      out += i->second->dot_product(*j.second);
      ++i;
    }
  } else {
    throw logic_error("Trying to copy something different in StorageIncore");
  }
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class StorageBlock<double>;
template class StorageBlock<complex<double>>;
template class StorageIncore<double>;
template class StorageIncore<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

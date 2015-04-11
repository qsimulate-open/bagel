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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename DataType>
Storage_Incore<DataType>::Storage_Incore(const map<size_t, size_t>& size, bool init) : Storage_base<StorageBlock<DataType>>(size, init) {
  static_assert(is_same<DataType, double>::value or is_same<DataType, complex<double>>::value, "illegal Type in Storage_Incore");
}


template<typename DataType>
unique_ptr<DataType[]> Storage_Incore<DataType>::get_block(const size_t& key) const {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::get_block(const size_t&)");
  return hash->second->get_block();
}


template<typename DataType>
unique_ptr<DataType[]> Storage_Incore<DataType>::move_block(const size_t& key) {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::move_block(const size_t&)");
  return hash->second->move_block();
}


template<typename DataType>
void Storage_Incore<DataType>::put_block(const size_t& key, unique_ptr<DataType[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");
  hash->second->put_block(move(dat));
}


template<typename DataType>
void Storage_Incore<DataType>::add_block(const size_t& key, const unique_ptr<DataType[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");
  hash->second->add_block(dat);
}


template<typename DataType>
void Storage_Incore<DataType>::zero() {
  for (auto& i : hashtable_)
    i.second->zero();
}


template<typename DataType>
void Storage_Incore<DataType>::scale(const DataType& a) {
  for (auto& i : hashtable_)
    i.second->scale(a);
}


template<typename DataType>
Storage_Incore<DataType>& Storage_Incore<DataType>::operator=(const Storage_Incore<DataType>& o) {
  if (hashtable_.size() == o.hashtable_.size()) {
    auto i = hashtable_.begin();
    for (auto& j : o.hashtable_) {
      *i->second = *j.second;
      ++i;
    }
  } else {
    throw logic_error("Trying to copy something different in Storage_Incore");
  }
  return *this;
}


template<typename DataType>
void Storage_Incore<DataType>::ax_plus_y(const DataType& a, const Storage_Incore<DataType>& o) {
  if (hashtable_.size() == o.hashtable_.size()) {
    auto i = hashtable_.begin();
    for (auto& j : o.hashtable_) {
      i->second->ax_plus_y(a, *j.second);
      ++i;
    }
  } else {
    throw logic_error("Trying to copy something different in Storage_Incore");
  }
}


template<typename DataType>
DataType Storage_Incore<DataType>::dot_product(const Storage_Incore<DataType>& o) const {
  DataType out = 0.0;
  if (hashtable_.size() == o.hashtable_.size()) {
    auto i = hashtable_.begin();
    for (auto& j : o.hashtable_) {
      out += i->second->dot_product(*j.second);
      ++i;
    }
  } else {
    throw logic_error("Trying to copy something different in Storage_Incore");
  }
  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class StorageBlock<double>;
template class StorageBlock<complex<double>>;
template class Storage_Incore<double>;
template class Storage_Incore<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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


StorageBlock::StorageBlock(const size_t size, const bool init) : size_(size), initialized_(init) {
  if (init) {
    data_ = unique_ptr<double[]>(new double[size_]);
    zero();
  }
}


void StorageBlock::zero() {
  if (initialized_)
    fill_n(data(), size_, 0.0);
}


StorageBlock& StorageBlock::operator=(const StorageBlock& o) {
  if (o.initialized_ && !initialized_) {
    data_ = unique_ptr<double[]>(new double[size_]);
    initialized_ = true;
  }
  if (o.initialized_)
    copy_n(o.data(), size_, data());
  return *this;
}


void StorageBlock::put_block(unique_ptr<double[]>&& o) {
  assert(!initialized_);
  initialized_ = true;
  data_ = move(o);
}


void StorageBlock::add_block(const std::unique_ptr<double[]>& o) {
  assert(initialized_);
  blas::ax_plus_y_n(1.0, o.get(), size_, data());
}


unique_ptr<double[]> StorageBlock::get_block() const {
  assert(initialized_);
  unique_ptr<double[]> out(new double[size_]);
  copy_n(data_.get(), size_, out.get());
  return move(out);
}


unique_ptr<double[]> StorageBlock::move_block() {
  if (!initialized_) {
    initialized_ = true;
    data_ = unique_ptr<double[]>(new double[size_]);
    zero();
  }
  initialized_ = false;
  return move(data_);
}


double StorageBlock::dot_product(const StorageBlock& o) const {
  assert(size_ == o.size_ && !(initialized_ ^ o.initialized_));
  return initialized_ ? blas::dot_product(data(), size_, o.data()) : 0.0;
}


void StorageBlock::ax_plus_y(const double a, const StorageBlock& o) {
  assert(size_ == o.size_ && !(initialized_ ^ o.initialized_));
  if (initialized_)
    blas::ax_plus_y_n(a, o.data(), size_, data());
}


void StorageBlock::scale(const double a) {
  if (initialized_)
    blas::scale_n(a, data(), size_);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Storage_Incore::Storage_Incore(const map<size_t, size_t>& size, bool init) : Storage_base<StorageBlock>(size, init) {
}


unique_ptr<double[]> Storage_Incore::get_block(const size_t& key) const {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::get_block(const size_t&)");
  return hash->second->get_block();
}


unique_ptr<double[]> Storage_Incore::move_block(const size_t& key) {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::move_block(const size_t&)");
  return hash->second->move_block();
}


void Storage_Incore::put_block(const size_t& key, unique_ptr<double[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");
  hash->second->put_block(move(dat));
}


void Storage_Incore::add_block(const size_t& key, const unique_ptr<double[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");
  hash->second->add_block(dat);
}


void Storage_Incore::zero() {
  for (auto& i : hashtable_)
    i.second->zero();
}


void Storage_Incore::scale(const double a) {
  for (auto& i : hashtable_)
    i.second->scale(a);
}


Storage_Incore& Storage_Incore::operator=(const Storage_Incore& o) {
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


void Storage_Incore::ax_plus_y(const double a, const Storage_Incore& o) {
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


double Storage_Incore::dot_product(const Storage_Incore& o) const {
  double out = 0.0;
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



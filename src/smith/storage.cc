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
#include <src/smith/storage.h>
#include <stdexcept>
#include <iostream>
#include <algorithm>

using namespace bagel::SMITH;
using namespace std;


Storage_Incore::Storage_Incore(const map<size_t, size_t>& size, bool init) : Storage_base(size, init) {
  if (init) {
    for (auto& i : size) {
      unique_ptr<double[]> tmp(new double[i.second]);
      fill(tmp.get(), tmp.get()+i.second, 0.0);
      data_.push_back(move(tmp));
    }
  } else {
    // if not init, we make a dummy tensor with size 1
    for (auto& i : size) {
      unique_ptr<double[]> tmp(new double[1]);
      data_.push_back(move(tmp));
    }
  }
}

void Storage_Incore::initialize() {
  for (auto& i : this->hashtable_) {
    unique_ptr<double[]> tmp(new double[i.second.second]);
    data_.push_back(move(tmp));
  }
}


unique_ptr<double[]> Storage_Incore::get_block(const size_t& key) const {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::get_block(const size_t&)");

  // then create a memory
  const size_t blocksize = hash->second.second;
  const size_t blocknum  = hash->second.first;
  unique_ptr<double[]> buf(new double[blocksize]);

  assert(initialized_[blocknum]);

  // then copy...
  copy_n(data_.at(blocknum).get(), blocksize, buf.get());

  return move(buf);
}


unique_ptr<double[]> Storage_Incore::move_block(const size_t& key) {
  // first find a key
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::move_block(const size_t&)");
  const size_t blocknum  = hash->second.first;

  if (!initialized(blocknum)) {
    unique_ptr<double[]> tmp(new double[hash->second.second]);
    fill_n(tmp.get(), hash->second.second, 0.0);
    data_[blocknum] = move(tmp);
    initialized_[blocknum] = true;
  }

  assert(initialized_[blocknum]);
  return move(data_.at(blocknum));
}


void Storage_Incore::put_block(const size_t& key, unique_ptr<double[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");
  const size_t blocknum  = hash->second.first;
  data_[blocknum] = move(dat);
  initialized_[blocknum] = true;
}


void Storage_Incore::add_block(const size_t& key, const unique_ptr<double[]>& dat) {
  auto hash = hashtable_.find(key);
  if (hash == hashtable_.end())
    throw logic_error("a key was not found in Storage::put_block(const size_t&)");

  const size_t blocksize = hash->second.second;
  const size_t blocknum  = hash->second.first;

  daxpy_(blocksize, 1.0, dat, 1, data_.at(blocknum), 1);
}


void Storage_Incore::zero() {
  for (auto& i : hashtable_) {
    const size_t bn = i.second.first;
    const size_t ln = i.second.second;
    fill(data_[bn].get(), data_[bn].get()+ln, 0.0);
  }
}

void Storage_Incore::scale(const double a) {
  for (auto& i : hashtable_) {
    const size_t bn = i.second.first;
    const size_t ln = i.second.second;
    dscal_(ln, a, data_[bn], 1);
  }
}

Storage_Incore& Storage_Incore::operator=(const Storage_Incore& o) {
  if (data_.size() == o.data_.size()) {
    auto i = data_.begin();
    auto k = hashtable_.begin();
    auto l = o.hashtable_.begin();
    for (auto j = o.data_.begin(); j != o.data_.end(); ++j, ++i, ++k, ++l) {
      if (k->second != l->second || k->first != l->first)
        throw logic_error("Trying to copy something different in Storage_Incore");
      copy_n(j->get(), k->second.second, i->get());
    }
  } else {
    throw logic_error("Trying to copy something different in Storage_Incore");
  }
  return *this;
}


void Storage_Incore::ax_plus_y(const double a, const Storage_Incore& o) {
  if (data_.size() == o.data_.size()) {
    auto i = data_.begin();
    auto k = hashtable_.begin();
    auto l = o.hashtable_.begin();
    for (auto j = o.data_.begin(); j != o.data_.end(); ++j, ++i, ++k, ++l) {
      if (k->second != l->second || k->first != l->first)
        throw logic_error("Trying to copy something different in Storage_Incore");
      daxpy_(k->second.second, a, *j, 1, *i, 1);
    }
  } else {
    throw logic_error("Trying to copy something different in Storage_Incore");
  }
}


double Storage_Incore::dot_product(const Storage_Incore& o) const {
  double out = 0.0;
  if (data_.size() == o.data_.size()) {
    auto i = data_.begin();
    auto k = hashtable_.begin();
    auto l = o.hashtable_.begin();
    for (auto j = o.data_.begin(); j != o.data_.end(); ++j, ++i, ++k, ++l) {
      if (k->second != l->second || k->first != l->first)
        throw logic_error("Trying to copy something different in Storage_Incore");
      out += ddot_(k->second.second, *j, 1, *i, 1);
    }
  } else {
    throw logic_error("Trying to copy something different in Storage_Incore");
  }
  return out;
}



//
// BAGEL - Parallel electron correlation program.
// Filename: storagekramers.cc
// Copyright (C) 2015 Toru Shiozaki
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

#include <src/smith/storagekramers.h>

using namespace bagel::SMITH;
using namespace std;


template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block() const {
  return get_block_();
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block(const Index& i0) const {
  return get_block_(i0);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block(const Index& i0, const Index& i1) const {
  return get_block_(i0, i1);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2) const {
  return get_block_(i0, i1, i2);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3) const {
  return get_block_(i0, i1, i2, i3);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4) const {
  return get_block_(i0, i1, i2, i3, i4);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5) const {
  return get_block_(i0, i1, i2, i3, i4, i5);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5, const Index& i6) const {
  return get_block_(i0, i1, i2, i3, i4, i5, i6);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                           const Index& i4, const Index& i5, const Index& i6, const Index& i7) const {
  return get_block_(i0, i1, i2, i3, i4, i5, i6, i7);
}


template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block() {
  return move_block_();
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block(const Index& i0) {
  return move_block_(i0);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block(const Index& i0, const Index& i1) {
  return move_block_(i0, i1);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2) {
  return move_block_(i0, i1, i2);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  return move_block_(i0, i1, i2, i3);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                            const Index& i4) {
  return move_block_(i0, i1, i2, i3, i4);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                            const Index& i4, const Index& i5) {
  return move_block_(i0, i1, i2, i3, i4, i5);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                            const Index& i4, const Index& i5, const Index& i6) {
  return move_block_(i0, i1, i2, i3, i4, i5, i6);
}

template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::move_block(const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                            const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  return move_block_(i0, i1, i2, i3, i4, i5, i6, i7);
}


template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat) {
  put_block_(dat);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0) {
  put_block_(dat, i0);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) {
  put_block_(dat, i0, i1);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) {
  put_block_(dat, i0, i1, i2);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  put_block_(dat, i0, i1, i2, i3);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                      const Index& i4) {
  put_block_(dat, i0, i1, i2, i3, i4);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                      const Index& i4, const Index& i5) {
  put_block_(dat, i0, i1, i2, i3, i4, i5);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                      const Index& i4, const Index& i5, const Index& i6) {
  put_block_(dat, i0, i1, i2, i3, i4, i5, i6);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                      const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  put_block_(dat, i0, i1, i2, i3, i4, i5, i6, i7);
}


template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat) {
  add_block_(dat);
}

template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0) {
  add_block_(dat, i0);
}

template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) {
  add_block_(dat, i0, i1);
}

template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) {
  add_block_(dat, i0, i1, i2);
}

template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  add_block_(dat, i0, i1, i2, i3);
}

template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                            const Index& i4) {
  add_block_(dat, i0, i1, i2, i3, i4);
}

template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                            const Index& i4, const Index& i5) {
  add_block_(dat, i0, i1, i2, i3, i4, i5);
}

template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                            const Index& i4, const Index& i5, const Index& i6) {
  add_block_(dat, i0, i1, i2, i3, i4, i5, i6);
}

template<typename DataType>
void StorageKramers<DataType>::add_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                            const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  add_block_(dat, i0, i1, i2, i3, i4, i5, i6, i7);
}

// explicit instantiation
template class StorageKramers<double>;
template class StorageKramers<complex<double>>;

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: storagekramers.cc
// Copyright (C) 2015 Toru Shiozaki
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

#include <src/smith/storagekramers.h>

using namespace bagel::SMITH;
using namespace std;


template<typename DataType>
unique_ptr<DataType[]> StorageKramers<DataType>::get_block() const {
  return RMAWindow<DataType>::rma_get(generate_hash_key());
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
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat) {
  put_block_(dat);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0) {
  put_block_(dat, i0);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1) {
  put_block_(dat, i0, i1);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2) {
  put_block_(dat, i0, i1, i2);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3) {
  put_block_(dat, i0, i1, i2, i3);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                            const Index& i4) {
  put_block_(dat, i0, i1, i2, i3, i4);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                            const Index& i4, const Index& i5) {
  put_block_(dat, i0, i1, i2, i3, i4, i5);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                            const Index& i4, const Index& i5, const Index& i6) {
  put_block_(dat, i0, i1, i2, i3, i4, i5, i6);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, const Index& i0, const Index& i1, const Index& i2, const Index& i3,
                                                                            const Index& i4, const Index& i5, const Index& i6, const Index& i7) {
  put_block_(dat, i0, i1, i2, i3, i4, i5, i6, i7);
}

template<typename DataType>
void StorageKramers<DataType>::put_block(const unique_ptr<DataType[]>& dat, vector<Index> indices) {
  put_block_(dat, indices);
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
template class bagel::SMITH::StorageKramers<double>;
template class bagel::SMITH::StorageKramers<complex<double>>;

BOOST_CLASS_EXPORT_IMPLEMENT(StorageKramers<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(StorageKramers<complex<double>>)

#endif

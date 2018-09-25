//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: tensor.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/smith/tensor.h>
#include <src/smith/storagekramers.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

template <typename DataType>
Tensor_<DataType>::Tensor_(vector<IndexRange> in, const bool kramers, const unordered_set<size_t> sparse, const bool alloc)
  : range_(in), rank_(in.size()), sparse_(sparse), initialized_(false), allocated_(alloc) {

  // make block list
  // First make sure the tensor is not empty
  if (!in.empty() && !any_of(in.begin(), in.end(), [](IndexRange i){return (i.range().size() == 0);})) {
    vector<vector<Index>> index = LoopGenerator::gen(in);

    // first compute hashtags and length
    map<size_t, size_t> hashmap;
    for (auto& i : index) {
      size_t size = 1lu;
      vector<size_t> h;
      for (auto& j : i) {
        size *= j.size();
        h.push_back(j.key());
      }
      auto key = generate_hash_key(h);
      if (sparse.empty() || sparse.count(key))
        hashmap.emplace(key, size);
    }

    if (!kramers)
      data_ = make_shared<Storage<DataType>>(hashmap, alloc);
    else
      data_ = make_shared<StorageKramers<DataType>>(hashmap, alloc);
  } else {
    rank_ = 0;
    map<size_t, size_t> hashmap {{generate_hash_key(), 1lu}};
    data_ = make_shared<Storage<DataType>>(hashmap, alloc);
  }
}


template <typename DataType>
void Tensor_<DataType>::allocate() {
  assert(!allocated_);
  data_->initialize();
  allocated_ = true;
}


template <typename DataType>
size_t Tensor_<DataType>::size_alloc() const {
  return data_->size_alloc();
}


template <typename DataType>
Tensor_<DataType>& Tensor_<DataType>::operator=(const Tensor_<DataType>& o) {
  *data_ = *(o.data_);
  allocated_ = true;
  return *this;
}


template <typename DataType>
shared_ptr<Tensor_<DataType>> Tensor_<DataType>::clone() const {
  return make_shared<Tensor_<DataType>>(range_, false, sparse_, true);
}


template <typename DataType>
shared_ptr<Tensor_<DataType>> Tensor_<DataType>::copy() const {
  shared_ptr<Tensor_<DataType>> out = clone();
  *out = *this;
  return out;
}


template <typename DataType>
vector<DataType> Tensor_<DataType>::diag() const {
  if (rank_ != 2 || range_.at(0) != range_.at(1))
    throw logic_error("Tensor_<DataType>::diag can be called only with a square tensor of rank 2");
  const size_t size = range_.at(0).back().offset() + range_.at(0).back().size();
  vector<DataType> buf(size);
  for (auto& i : range_.at(0)) {
    unique_ptr<DataType[]> data0 = get_block(i, i);
    for (int j = 0; j != i.size(); ++j) {
      buf[j+i.offset()] = data0[j+j*i.size()];
    }
  }
  return buf;
}


template <typename DataType>
double Tensor_<DataType>::orthog(const list<shared_ptr<const Tensor_<DataType>>> o) {
  for (auto& it : o)
    ax_plus_y(-detail::conj(this->dot_product(it)), it);
  const double n = norm();
  scale(1.0/n);
  return n;
}


template <typename DataType>
shared_ptr<typename Tensor_<DataType>::MatType> Tensor_<DataType>::matrix() const {
  vector<IndexRange> o = indexrange();
  assert(o.size() == 2);

  const int off0 = o[0].front().offset();
  const int off1 = o[1].front().offset();
  auto out = make_shared<MatType>(o[0].size(), o[1].size());
  for (auto& i1 : o[1].range())
    for (auto& i0 : o[0].range())
      out->copy_block(i0.offset()-off0, i1.offset()-off1, i0.size(), i1.size(), get_block(i0, i1).get());

  return out;
}


template <typename DataType>
shared_ptr<typename Tensor_<DataType>::VecType> Tensor_<DataType>::vectorb() const {
  vector<IndexRange> o = indexrange();
  assert(o.size() == 1);

  auto out = make_shared<VecType>(o[0].size());
  for (auto& i0 : o[0].range())
    copy_n(get_block(i0).get(), i0.size(), out->data()+i0.offset());

  return out;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class SMITH::Tensor_<double>;
template class SMITH::Tensor_<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_CLASS_EXPORT_IMPLEMENT(Tensor_<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(Tensor_<complex<double>>)

#endif

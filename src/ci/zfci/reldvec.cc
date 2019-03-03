//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reldvec.cc
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

#include <numeric>
#include <src/ci/zfci/reldvec.h>

using namespace std;
using namespace bagel;

template<typename DataType>
RelDvector<DataType>::RelDvector(shared_ptr<const Space_base> space, const size_t ij) : space_(space) {
  for (auto& isp : space->detmap())
    dvecs_.emplace(make_pair(isp.second->nelea(), isp.second->neleb()),
                   make_shared<Dvector<DataType>>(isp.second, ij));
}


template<typename DataType>
RelDvector<DataType>::RelDvector(const RelDvector<DataType>& o) : space_(o.space_) {
  for (auto& i : o.dvecs_)
    dvecs_.emplace(i.first, i.second->copy());
}


template<typename DataType>
RelDvector<DataType>::RelDvector(const vector<shared_ptr<RelDvector<DataType>>>& o) : space_(o.front()->space_) {
  for (auto& isp : space_->detmap())
    dvecs_.emplace(make_pair(isp.second->nelea(), isp.second->neleb()), make_shared<Dvector<DataType>>(isp.second, o.size()));
  int j = 0;
  for (auto& i : o)
    set_data(j++, i);
}


template<typename DataType>
shared_ptr<RelDvector<DataType>> RelDvector<DataType>::clone() const {
  return make_shared<RelDvector<DataType>>(space_, dvecs_.begin()->second->ij());
}


template<typename DataType>
shared_ptr<RelDvector<DataType>> RelDvector<DataType>::copy() const {
  return make_shared<RelDvector<DataType>>(*this);
}


template<typename DataType>
shared_ptr<RelDvector<DataType>> RelDvector<DataType>::extract_state(const vector<int> input) const {
  MapType newdvec;
  for (auto& i : dvecs_)
    newdvec.emplace(i.first, i.second->extract_state(input));
  return make_shared<RelDvector<DataType>>(newdvec, space_);
}


template<typename DataType>
void RelDvector<DataType>::set_data(const int istate, shared_ptr<const RelDvector<DataType>> o) {
  assert(space_ == o->space_ || o->dvecs_.begin()->second->ij() == 1);
  auto j = o->dvecs_.begin();
  for (auto& i : dvecs_) {
    *i.second->data(istate) = *j->second->data(0);
    ++j;
  }
}


template<typename DataType>
void RelDvector<DataType>::zero() {
  for_each(dvecs_.begin(), dvecs_.end(), [](EleType i) { i.second->zero(); });
}


template<typename DataType>
size_t RelDvector<DataType>::size() const {
  return accumulate(dvecs_.begin(), dvecs_.end(), 0ull, [](size_t i, const EleType& o) { return i+o.second->size(); });
}


template<typename DataType>
DataType RelDvector<DataType>::dot_product(const RelDvector<DataType>& o) const {
  return inner_product(dvecs_.begin(), dvecs_.end(), o.dvecs_.begin(), DataType(0.0), plus<DataType>(),
                            [](const EleType& i, const EleType& j) { return i.second->dot_product(*j.second); });
}


template<typename DataType>
void RelDvector<DataType>::ax_plus_y(const DataType a, const RelDvector<DataType>& o) {
  auto iter = o.dvecs_.begin();
  for (auto& i : dvecs_) {
    assert(i.first == iter->first);
    i.second->ax_plus_y(a, *iter->second);
    ++iter;
  }
}


template<typename DataType>
vector<shared_ptr<const RelDvector<DataType>>> RelDvector<DataType>::split(const int nstart, const int nend) const {
  vector<shared_ptr<const RelDvector<DataType>>> out;
  for (int i = nstart; i != nend; ++i) {
    MapType tmp;
    // copy construct each of them
    for (auto& j : dvecs_) {
      vector<shared_ptr<Civector<DataType>>> tmp1 { make_shared<Civector<DataType>>(*j.second->data(i)) };
      tmp.emplace(j.first, make_shared<Dvector<DataType>>(tmp1));
    }
    out.push_back(make_shared<RelDvector<DataType>>(tmp, space_));
  }
  return out;
}


template<typename DataType>
vector<shared_ptr<const RelDvector<DataType>>> RelDvector<DataType>::dvec(const vector<int>& conv) const {
  vector<shared_ptr<const RelDvector<DataType>>> sp = split();
  vector<shared_ptr<const RelDvector<DataType>>> out;
  auto c = conv.begin();
  for (auto& i : sp)
    if (*c++ == 0) out.push_back(i);
    else out.push_back(nullptr);
  return out;
}


template<typename DataType>
void RelDvector<DataType>::scale(const DataType& a) {
  for_each(dvecs_.begin(), dvecs_.end(), [&a](EleType i) { i.second->scale(a); });
}


template<typename DataType>
double RelDvector<DataType>::orthog(list<shared_ptr<const RelDvector<DataType>>> c) {
  for (auto& iter : c)
    project_out(iter);
  return normalize();
}


template<typename DataType>
double RelDvector<DataType>::orthog(shared_ptr<const RelDvector<DataType>> o) {
  return orthog(list<shared_ptr<const RelDvector<DataType>>>{o});
}


template<typename DataType>
double RelDvector<DataType>::normalize() {
  const double norm = this->norm();
  const double scal = (norm*norm<1.0e-60 ? 0.0 : 1.0/norm);
  scale(DataType(scal));
  return norm;
}


template<typename DataType>
void RelDvector<DataType>::print(double thresh) const {
  for (auto& i : dvecs_)
    i.second->print(thresh);
}


template<typename DataType>
void RelDvector<DataType>::synchronize() {
#ifdef HAVE_MPI_H
  for (auto& i : dvecs_)
    i.second->synchronize();
#endif
}


template class bagel::RelDvector<double>;
template class bagel::RelDvector<complex<double>>;

BOOST_CLASS_EXPORT_IMPLEMENT(RelDvector<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(RelDvector<complex<double>>)

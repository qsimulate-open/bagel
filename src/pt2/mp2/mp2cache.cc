//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mp2cache.h
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

#include <src/pt2/mp2/mp2cache.h>

using namespace std;
using namespace bagel;

template<>
MP2Tag<double> MP2Cache_<double>::create_cache_and_request_recv(const int i, const int origin) {
  cache_[i] = make_shared<Matrix>(naux_, nvirt_, true);
  return MP2Tag<double>{mpi__->request_recv(cache_[i]->data(), cache_[i]->size(), origin, myrank_*nocc_+i)};
}

template<>
MP2Tag<complex<double>> MP2Cache_<complex<double>>::create_cache_and_request_recv(const int i, const int origin) {
  cache_[i] = {make_shared<Matrix>(naux_, nvirt_, true), make_shared<Matrix>(naux_, nvirt_, true)};
  return MP2Tag<complex<double>>{mpi__->request_recv(cache_[i].first->data(),  cache_[i].first->size(),  origin, myrank_*nocc_*2+i*2),
                                 mpi__->request_recv(cache_[i].second->data(), cache_[i].second->size(), origin, myrank_*nocc_*2+i*2+1)};
}

template<>
MP2Tag<double> MP2Cache_<double>::request_send(const int i, const int dest) {
  return MP2Tag<double>{mpi__->request_send(fullt_->data() + (i*nvirt_-fullt_->bstart())*naux_, nvirt_*naux_, dest, dest*nocc_+i)};
}

template<>
MP2Tag<complex<double>> MP2Cache_<complex<double>>::request_send(const int i, const int dest) {
  return MP2Tag<complex<double>>{mpi__->request_send(fullt_->data(0) + (i*nvirt_-fullt_->bstart())*naux_, nvirt_*naux_, dest, dest*nocc_*2+i*2),
                                 mpi__->request_send(fullt_->data(1) + (i*nvirt_-fullt_->bstart())*naux_, nvirt_*naux_, dest, dest*nocc_*2+i*2+1)};
}

template<>
shared_ptr<const Matrix> MP2Cache_<double>::operator()(const int i) const {
  return cache_.at(i);
}

template<>
shared_ptr<const ZMatrix> MP2Cache_<complex<double>>::operator()(const int i) const {
  return make_shared<ZMatrix>(*cache_.at(i).first, *cache_.at(i).second);
}

/// specialization ///
template class bagel::MP2Cache_<double>;
template class bagel::MP2Cache_<complex<double>>;

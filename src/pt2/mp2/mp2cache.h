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

#ifndef __SRC_PT2_MP2_MP2CACHE_H
#define __SRC_PT2_MP2_MP2CACHE_H

#include <set>
#include <src/df/dfdistt.h>
#include <src/df/reldffullt.h>

namespace bagel {

template <typename DataType>
class MP2Tag {
  protected:
    typename std::conditional<std::is_same<DataType,double>::value,std::array<int,1>,std::array<int,2>>::type tag_;

  public:
    MP2Tag() { std::fill(tag_.begin(), tag_.end(), -1); }
    MP2Tag(std::initializer_list<int> o) {
      assert(o.size() == tag_.size());
      auto j = tag_.begin();
      for (auto& i : o)
        *j++ = i;
    }
    MP2Tag(const MP2Tag<DataType>&) = default;
    MP2Tag(MP2Tag<DataType>&&) = default;

    MP2Tag<DataType>& operator=(const MP2Tag<DataType>& o) = default; 
    MP2Tag<DataType>& operator=(MP2Tag<DataType>&&) = default;

    bool operator<(const MP2Tag<DataType>& o) const {
      auto j = tag_.begin();
      for (auto i : o.tag_) {
        if (*j != i) return *j < i;
        ++j;
      }
      return false;
    }

    bool operator==(const MP2Tag<DataType>& o) const {
      return !(*this < o) && !(o < *this);
    }
    bool invalid() const { return tag_.at(0) == -1; }
    void wait() const {
      if (!invalid())
        for (auto& i : tag_)
          mpi__->wait(i);
    }
    bool test() const {
      bool out = true;
      for (auto& i : tag_) out &= mpi__->test(i);
      return out;
    }
};

template <typename DataType>
class MP2Cache_ {
  protected:
    using CacheType = typename std::conditional<std::is_same<DataType,double>::value, std::shared_ptr<Matrix>,
                                                                                      std::pair<std::shared_ptr<Matrix>,std::shared_ptr<Matrix>>>::type;
    using DFType    = typename std::conditional<std::is_same<DataType,double>::value, DFDistT, RelDFFullT>::type;
    using MatType   = typename std::conditional<std::is_same<DataType,double>::value, Matrix, ZMatrix>::type;

  protected:
    const size_t naux_;
    const size_t nocc_;
    const size_t nvirt_;

    const std::shared_ptr<const DFType> fullt_; // aux, virt, occ
    std::vector<std::vector<std::tuple<int,int,MP2Tag<DataType>,MP2Tag<DataType>>>> tasks_;

    // the data is stored in a map
    std::map<int, CacheType> cache_;
    // pair of node and set of integers
    std::vector<std::set<int>> cachetable_;
    std::vector<MP2Tag<DataType>> sendreqs_;

    int myrank_;
    int nloop_;

    MP2Tag<DataType> create_cache_and_request_recv(const int i, const int origin);
    MP2Tag<DataType> request_one(const int i, const int rank) {
      if (i < 0 || i >= nocc_) return MP2Tag<DataType>();
      cachetable_[rank].insert(i);
      MP2Tag<DataType> tag;
      if (cache_.find(i) == cache_.end() && myrank_ == rank) {
        const int origin = fullt_->locate(0, i*nvirt_);
        if (origin == myrank_)
          cache_[i] = fullt_->get_slice(i*nvirt_, (i+1)*nvirt_).front();
        else
          tag = create_cache_and_request_recv(i, origin);
      }
      return tag;
    }

    MP2Tag<DataType> request_send(const int i, const int dest);
    MP2Tag<DataType> send_one(const int i, const int dest) {
      // see if "i" is cached at dest
      if (i < 0 || i >= nocc_ || cachetable_[dest].count(i) || fullt_->locate(0, i*nvirt_) != myrank_)
        return MP2Tag<DataType>();
      return request_send(i, dest);
    }

  public:
    MP2Cache_(const int naux, const int nocc, const int nvirt, std::shared_ptr<const DFType> fullt,
              const std::vector<std::vector<std::tuple<int,int,MP2Tag<DataType>,MP2Tag<DataType>>>>& tasks
                  = std::vector<std::vector<std::tuple<int,int,MP2Tag<DataType>,MP2Tag<DataType>>>>())
     : naux_(naux), nocc_(nocc), nvirt_(nvirt), fullt_(fullt), tasks_(tasks), cachetable_(mpi__->size()), myrank_(mpi__->rank()) {

      assert(naux_ == fullt->naux());

      // if not specified make a list for static distribution of ij
      if (tasks_.empty()) {
        tasks_.resize(mpi__->size());
        int nmax = 0;
        StaticDist ijdist(nocc*(nocc+1)/2, mpi__->size());
        for (int inode = 0; inode != mpi__->size(); ++inode) {
          for (int i = 0, cnt = 0; i < nocc; ++i)
            for (int j = i; j < nocc; ++j, ++cnt)
              if (cnt >= ijdist.start(inode) && cnt < ijdist.start(inode) + ijdist.size(inode))
                tasks_[inode].push_back(std::make_tuple(j, i, MP2Tag<DataType>(),MP2Tag<DataType>()));
          if (tasks_[inode].size() > nmax) nmax = tasks_[inode].size();
        }
        for (auto& i : tasks_) {
          const int n = i.size();
          for (int j = 0; j != nmax-n; ++j) i.push_back(std::make_tuple(-1,-1,MP2Tag<DataType>(),MP2Tag<DataType>()));
        }
      }
      nloop_ = tasks_[0].size();
    }

    std::shared_ptr<const MatType> operator()(const int i) const;

    const std::tuple<int,int,MP2Tag<DataType>,MP2Tag<DataType>>& task(const int i) const { return tasks_[myrank_][i]; }
    const std::vector<std::vector<std::tuple<int,int,MP2Tag<DataType>,MP2Tag<DataType>>>>& tasks() const { return tasks_; }

    int nloop() const { return nloop_; }

    void block(const int nadd, const int ndrop) {
      assert(ndrop < nadd);
      if (ndrop >= 0) {
        for (int inode = 0; inode != mpi__->size(); ++inode) {
          const int id = std::get<0>(tasks_[inode][ndrop]);
          const int jd = std::get<1>(tasks_[inode][ndrop]);
          // if id and jd are no longer used in the cache, delete the element
          std::set<int> used;
          for (int i = ndrop+1; i <= nadd; ++i) {
            used.insert(std::get<0>(tasks_[inode][i]));
            used.insert(std::get<1>(tasks_[inode][i]));
          }
          if (id >= 0 && id < nocc_ && !used.count(id)) {
            if (inode == myrank_) cache_.erase(id);
            cachetable_[inode].erase(id);
          }
          if (jd >= 0 && jd < nocc_ && !used.count(jd)) {
            if (inode == myrank_) cache_.erase(jd);
            cachetable_[inode].erase(jd);
          }
        }
      }
      if (nadd < nloop_) {
        for (int inode = 0; inode != mpi__->size(); ++inode) {
          if (inode == myrank_) {
            // recieve data from other processes
            std::get<2>(tasks_[myrank_][nadd]) = request_one(std::get<0>(tasks_[myrank_][nadd]), myrank_); // receive requests
            std::get<3>(tasks_[myrank_][nadd]) = request_one(std::get<1>(tasks_[myrank_][nadd]), myrank_);
          } else {
            // send data to other processes
            const MP2Tag<DataType> i = send_one(std::get<0>(tasks_[inode][nadd]), inode); // send requests
            if (!i.invalid()) sendreqs_.push_back(i);
            request_one(std::get<0>(tasks_[inode][nadd]), inode); // update cachetable_
            const MP2Tag<DataType> j = send_one(std::get<1>(tasks_[inode][nadd]), inode); // send requests
            if (!j.invalid()) sendreqs_.push_back(j);
            request_one(std::get<1>(tasks_[inode][nadd]), inode);
          }
        }
      }
    }

    void data_wait(const int n) const {
      const MP2Tag<DataType> ti = std::get<2>(task(n));
      const MP2Tag<DataType> tj = std::get<3>(task(n));
      ti.wait();
      tj.wait();
    }

    void wait() const {
      for (auto& i : sendreqs_)
        i.wait();
    }

};

template<> MP2Tag<double> MP2Cache_<double>::create_cache_and_request_recv(const int i, const int origin);
template<> MP2Tag<double> MP2Cache_<double>::request_send(const int i, const int dest);
template<> std::shared_ptr<const Matrix> MP2Cache_<double>::operator()(const int i) const;

template<> MP2Tag<std::complex<double>> MP2Cache_<std::complex<double>>::create_cache_and_request_recv(const int i, const int origin);
template<> MP2Tag<std::complex<double>> MP2Cache_<std::complex<double>>::request_send(const int i, const int dest);
template<> std::shared_ptr<const ZMatrix> MP2Cache_<std::complex<double>>::operator()(const int i) const;

extern template class MP2Cache_<double>;
extern template class MP2Cache_<std::complex<double>>;

using MP2Cache = MP2Cache_<double>;
using RelMP2Cache = MP2Cache_<std::complex<double>>;

}

#endif

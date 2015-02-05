//
// BAGEL - Parallel electron correlation program.
// Filename: mp2cache.h
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

#ifndef __SRC_PT2_MP2_MP2CACHE_H
#define __SRC_PT2_MP2_MP2CACHE_H

#include <set>
#include <src/df/dfdistt.h>

namespace bagel {

class MP2Cache {
  protected:
    const size_t naux_;
    const size_t nocc_;
    const size_t nvirt_;

    const std::shared_ptr<const DFDistT> fullt_;
    std::vector<std::vector<std::tuple<int,int,int,int>>> tasks_;

    // the data is stored in a map
    std::map<int, std::shared_ptr<Matrix>> cache_;
    // pair of node and set of integers
    std::vector<std::set<int>> cachetable_;
    std::vector<int> sendreqs_;

    const int myrank_;
    const int nloop_;

  public:
    MP2Cache(const int naux, const int nocc, const int nvirt, std::shared_ptr<const DFDistT> fullt, const std::vector<std::vector<std::tuple<int,int,int,int>>>& tasks)
     : naux_(naux), nocc_(nocc), nvirt_(nvirt), fullt_(fullt), tasks_(tasks), cachetable_(mpi__->size()), myrank_(mpi__->rank()), nloop_(tasks_[0].size()) {
    }

    std::shared_ptr<const Matrix> operator()(const int i) const { return cache_.at(i); }

    int nloop() const { return nloop_; }

    const std::tuple<int,int,int,int>& task(const int i) const { return tasks_[myrank_][i]; }

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
          if (!used.count(id)) {
            if (inode == myrank_) cache_.erase(id);
            cachetable_[inode].erase(id);
          }
          if (!used.count(jd)) {
            if (inode == myrank_) cache_.erase(jd);
            cachetable_[inode].erase(jd);
          }
        }
      }
      if (nadd < nloop_) {
        // issue recv requests
        auto request_one_ = [&](const int i, const int rank) {
          if (i < 0 || i >= nocc_) return -1;
          cachetable_[rank].insert(i);
          int tag = -1;
          if (cache_.find(i) == cache_.end() && myrank_ == rank) {
            const int origin = fullt_->locate(0, i*nvirt_);
            if (origin == myrank_) {
              cache_[i] = fullt_->get_slice(i*nvirt_, (i+1)*nvirt_).front();
            } else {
              cache_[i] = std::make_shared<Matrix>(naux_, nvirt_, true);
              tag = mpi__->request_recv(cache_[i]->data(), cache_[i]->size(), origin, myrank_*nocc_+i);
            }
          }
          return tag;
        };

        // issue send requests
        auto send_one_ = [&](const int i, const int dest) {
          // see if "i" is cached at dest
          if (i < 0 || i >= nocc_ || cachetable_[dest].count(i) || fullt_->locate(0, i*nvirt_) != myrank_)
            return -1;
          return mpi__->request_send(fullt_->data() + (i*nvirt_-fullt_->bstart())*naux_, nvirt_*naux_, dest, dest*nocc_+i);
        };

        for (int inode = 0; inode != mpi__->size(); ++inode) {
          if (inode == myrank_) {
            // recieve data from other processes
            std::get<2>(tasks_[myrank_][nadd]) = request_one_(std::get<0>(tasks_[myrank_][nadd]), myrank_); // receive requests
            std::get<3>(tasks_[myrank_][nadd]) = request_one_(std::get<1>(tasks_[myrank_][nadd]), myrank_);
          } else {
            // send data to other processes
            const int i = send_one_(std::get<0>(tasks_[inode][nadd]), inode); // send requests
            if (i >= 0) sendreqs_.push_back(i);
            request_one_(std::get<0>(tasks_[inode][nadd]), inode); // update cachetable_
            const int j = send_one_(std::get<1>(tasks_[inode][nadd]), inode); // send requests
            if (j >= 0) sendreqs_.push_back(j);
            request_one_(std::get<1>(tasks_[inode][nadd]), inode);
          }
        }
      }
    }

    void wait() const {
      for (auto& i : sendreqs_)
        mpi__->wait(i);
    }

};

}

#endif

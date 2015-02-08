//
// BAGEL - Parallel electron correlation program.
// Filename: mp2accum.h
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

#ifndef __SRC_PT2_MP2_MP2ACCUM_H
#define __SRC_PT2_MP2_MP2ACCUM_H

#include <set>
#include <src/df/dfdistt.h>

namespace bagel {

// TODO maybe good to make a base of MP2Accum and MP2Cache
class MP2Accum {
  protected:
    const size_t naux_;
    const size_t nocc_;
    const size_t nvirt_;

    const std::shared_ptr<DFDistT> fullt_; // aux, virt, occ
    std::vector<std::vector<std::tuple<int,int,int,int>>> tasks_;

    std::map<int, std::shared_ptr<const Matrix>> sendreqs_; // mpi tag and buffer
    std::map<int, std::pair<std::shared_ptr<Matrix>,int>> recvreqs_; // mpi tag, buffer and offset

    const int myrank_;
    const int nloop_;

  public:
    MP2Accum(const int naux, const int nocc, const int nvirt, std::shared_ptr<DFDistT> fullt,
             const std::vector<std::vector<std::tuple<int,int,int,int>>>& tasks)
     : naux_(naux), nocc_(nocc), nvirt_(nvirt), fullt_(fullt), tasks_(tasks), myrank_(mpi__->rank()), nloop_(tasks_[0].size()) {
      assert(naux_ == fullt->naux() && nvirt_ == fullt->nindex1() && nocc_ == fullt->nindex2());
    }

    template <int N, class = typename std::enable_if<N == 0 or N == 1>::type>
    void accumulate(const int n, std::shared_ptr<const Matrix> send) {
      assert(send->ndim() == naux_ && send->mdim() == nvirt_);
      // first scan sendreqs and recvreqs and drop matrices that have been sent/received
      for (auto i = sendreqs_.begin(); i != sendreqs_.end(); ) {
        if (mpi__->test(i->first))
          sendreqs_.erase(i++);
        else
          ++i;
      }
      for (auto i = recvreqs_.begin(); i != recvreqs_.end(); ) {
        if (mpi__->test(i->first)) {
          blas::ax_plus_y_n(1.0, i->second.first->data(), i->second.first->size(), fullt_->data()+fullt_->offset(0, i->second.second*nvirt_));
          recvreqs_.erase(i++);
        } else {
          ++i;
        }
      }

      // tag should not overlap with MP2Cache
      auto tag = [&](const int dest, const int origin) { return n + nloop_*(dest + mpi__->size()*origin) + nocc_*mpi__->size(); };

      // issue send requests
      auto send_one_ = [&](const int i) {
        if (i < 0 || i >= nocc_) return;
        const int dest = fullt_->locate(0, i*nvirt_);
        if (dest == myrank_) {
          blas::ax_plus_y_n(1.0, send->data(), send->size(), fullt_->data()+fullt_->offset(0, i*nvirt_));
        } else {
          sendreqs_.emplace(mpi__->request_send(send->data(), send->size(), dest, tag(dest, myrank_)), send);
        }
      };

      // issue send requests
      auto recv_one_ = [&](const int i) {
        if (i < 0 || i >= nocc_ || fullt_->locate(0, i*nvirt_) != myrank_)
          return;
        auto buf = std::make_shared<Matrix>(naux_, nvirt_);
        const int origin = fullt_->locate(0, i*nvirt_);
        recvreqs_.emplace(mpi__->request_send(buf->data(), buf->size(), origin, tag(myrank_, origin)), std::make_pair(buf, i));
      };

      // send, or locally accumulate
      send_one_(std::get<N>(tasks_[myrank_][n]));
      for (int inode = 0; inode != mpi__->size(); ++inode) {
        // check if somebody is sending data to me
        if (inode != myrank_)
          recv_one_(std::get<N>(tasks_[inode][n]));
      }
    }

    void wait() const {
      for (auto& i : recvreqs_) {
        mpi__->wait(i.first);
        blas::ax_plus_y_n(1.0, i.second.first->data(), i.second.first->size(), fullt_->data()+fullt_->offset(0, i.second.second*nvirt_));
      }
    }

};

}

#endif

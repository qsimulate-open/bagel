//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mp2accum.h
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

#ifndef __SRC_PT2_MP2_MP2ACCUM_H
#define __SRC_PT2_MP2_MP2ACCUM_H

#include <set>
#include <src/df/dfdistt.h>

namespace bagel {

class MP2Accum {
  protected:
    using Tag = MP2Tag<double>;
    const size_t naux_;
    const size_t nocc_;
    const size_t nvirt_;

    const std::shared_ptr<DFDistT> fullt_; // aux, virt, occ
    std::vector<std::vector<std::tuple<int,int,Tag,Tag>>> tasks_;

    std::map<Tag, std::shared_ptr<const Matrix>> sendreqs_; // mpi tag and buffer
    std::map<Tag, std::pair<std::shared_ptr<Matrix>,int>> recvreqs_; // mpi tag, buffer and offset

    const int myrank_;
    const int nloop_;

  public:
    MP2Accum(const int naux, const int nocc, const int nvirt, std::shared_ptr<DFDistT> fullt,
             const std::vector<std::vector<std::tuple<int,int,Tag,Tag>>>& tasks)
     : naux_(naux), nocc_(nocc), nvirt_(nvirt), fullt_(fullt), tasks_(tasks), myrank_(mpi__->rank()), nloop_(tasks_[0].size()) {
      assert(naux_ == fullt->naux() && nvirt_ == fullt->nindex1() && nocc_ == fullt->nindex2());
    }

    template <int N, class = typename std::enable_if<N == 0 or N == 1>::type>
    void accumulate(const int n, std::shared_ptr<const Matrix> send) {
      // first scan sendreqs and recvreqs and drop matrices that have been sent/received
      for (auto i = sendreqs_.begin(); i != sendreqs_.end(); ) {
        if (i->first.test())
          sendreqs_.erase(i++);
        else
          ++i;
      }
      for (auto i = recvreqs_.begin(); i != recvreqs_.end(); ) {
        if (i->first.test()) {
          blas::ax_plus_y_n(1.0, i->second.first->data(), i->second.first->size(), fullt_->data()+fullt_->offset(0, i->second.second*nvirt_));
          recvreqs_.erase(i++);
        } else {
          ++i;
        }
      }

      // tag should not overlap with MP2Cache
      auto tag = [&](const int dest, const int origin) { return N + 2*(n + nloop_*(dest + mpi__->size()*origin)) + nocc_*mpi__->size(); };

      // issue send requests
      auto send_one_ = [&](const int i) {
        if (i < 0 || i >= nocc_) return;
        const int dest = fullt_->locate(0, i*nvirt_);
        if (dest == myrank_) {
          blas::ax_plus_y_n(1.0, send->data(), send->size(), fullt_->data()+fullt_->offset(0, i*nvirt_));
        } else {
          sendreqs_.emplace(Tag{mpi__->request_send(send->data(), send->size(), dest, tag(dest, myrank_))}, send);
        }
      };

      // issue recv requests
      auto recv_one_ = [&](const int i, const int origin) {
        if (i < 0 || i >= nocc_ || fullt_->locate(0, i*nvirt_) != myrank_)
          return;
        auto buf = std::make_shared<Matrix>(naux_, nvirt_);
        recvreqs_.emplace(Tag{mpi__->request_recv(buf->data(), buf->size(), origin, tag(myrank_, origin))}, std::make_pair(buf, i));
      };

      // send, or locally accumulate
      if (N == 0 || std::get<0>(tasks_[myrank_][n]) != std::get<1>(tasks_[myrank_][n]))
        send_one_(std::get<N>(tasks_[myrank_][n]));

      for (int inode = 0; inode != mpi__->size(); ++inode) {
        // check if somebody is sending data to me
        if (inode != myrank_ && (N == 0 || std::get<0>(tasks_[inode][n]) != std::get<1>(tasks_[inode][n])))
          recv_one_(std::get<N>(tasks_[inode][n]), inode);
      }
    }

    void wait() {
      for (auto& i : recvreqs_) {
        i.first.wait();
        blas::ax_plus_y_n(1.0, i.second.first->data(), i.second.first->size(), fullt_->data()+fullt_->offset(0, i.second.second*nvirt_));
      }
      for (auto& i : sendreqs_)
        i.first.wait();
      recvreqs_.clear();
      sendreqs_.clear();
    }

};

}

#endif

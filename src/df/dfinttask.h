//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dfinttask.h
// Copyright (C) 2012 Toru Shiozaki
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

#ifndef __SRC_DF_DFINTTASK_H
#define __SRC_DF_DFINTTASK_H

#include <src/df/dfblock.h>
#include <src/molecule/shell.h>

namespace bagel {

template <typename TBatch, int N>
class DFIntTask {
  protected:
    const std::array<std::shared_ptr<const Shell>,4> shell_;
    const std::array<int,3> offset_; // at most 3 elements
    std::array<std::shared_ptr<DFBlock>,N> dfblocks_;

    std::shared_ptr<TBatch> compute_batch(const std::array<std::shared_ptr<const Shell>,4>& input) const {
      auto eribatch = std::make_shared<TBatch>(input, 2.0);
      eribatch->compute();
      return eribatch;
    }

  public:
    DFIntTask(std::array<std::shared_ptr<const Shell>,4>&& a, std::array<int,3>&& b, std::array<std::shared_ptr<DFBlock>,N>& df)
     : shell_(a), offset_(b), dfblocks_(df) { };

    void compute() {
      std::shared_ptr<TBatch> p = compute_batch(shell_);

      // all slot in
      for (int i = 0; i != N; ++i) {
        assert(dfblocks_[i]->b1size() == dfblocks_[i]->b2size());
        const size_t nbin = dfblocks_[i]->b1size();
        const size_t naux = dfblocks_[i]->asize();
        const double* ppt = p->data(i);
        double* const data = dfblocks_[i]->data();
        for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {
          for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1, ppt += shell_[1]->nbasis()) {
            std::copy_n(ppt, shell_[1]->nbasis(), data+offset_[2]+naux*(j1+nbin*j0));
            if (N == 1)
              std::copy_n(ppt, shell_[1]->nbasis(), data+offset_[2]+naux*(j0+nbin*j1));
          }
        }
      }
    }
};


}

#endif

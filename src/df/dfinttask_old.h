//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dfintask_old.h
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

#ifndef __SRC_DF_DFINTTASK_OLD_H
#define __SRC_DF_DFINTTASK_OLD_H

#include <src/integral/rys/rysintegral.h>

namespace bagel {

// T will be a DFDist
template <typename T>
class DFIntTask_OLD {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<int,2> offset_; // at most 3 elements
    int rank_;
    T* df_;

  public:
    DFIntTask_OLD(std::array<std::shared_ptr<const Shell>,4>&& a, std::array<int,2>&& b, T* df)
     : shell_(a), offset_(b), rank_(offset_.size()), df_(df) { }

    void compute() {

      std::pair<const double*, std::shared_ptr<RysIntegral<double, Int_t::Standard>>> p = df_->compute_batch(shell_);
      const double* ppt = p.first;

      const size_t naux = df_->naux();
      // all slot in
      if (rank_ == 2) {
        double* const data = df_->data2_->data();
        for (int j0 = offset_[0]; j0 != offset_[0] + shell_[2]->nbasis(); ++j0)
          for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++ppt)
            data[j1+j0*naux] = data[j0+j1*naux] = *ppt;
      } else {
        assert(false);
      }
    };
};

}

#endif

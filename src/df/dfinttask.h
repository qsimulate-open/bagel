//
// BAGEL - Parallel electron correlation program.
// Filename: dfinttask.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __SRC_DF_DFINTTASK_H
#define __SRC_DF_DFINTTASK_H

#include <src/scf/shell.h>
#include <src/rysint/rysint.h>

namespace bagel {

class DFIntTask {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<int,3> offset_; // at most 3 elements
    int rank_;
    DFBlock* dfblock_;

  public:
    DFIntTask(std::array<std::shared_ptr<const Shell>,4>&& a, std::vector<int>&& b, DFBlock* df) : shell_(a), rank_(b.size()), dfblock_(df) {
      int j = 0;
      for (auto i = b.begin(); i != b.end(); ++i, ++j) offset_[j] = *i;
    }
    DFIntTask() {};

    void compute() {
      std::pair<const double*, std::shared_ptr<Integral> > p = dfblock_->compute_batch(shell_);
      const double* ppt = p.first;

      assert(dfblock_->b1size_ == dfblock_->b2size_);
      const size_t nbin = dfblock_->b1size_;
      const size_t naux = dfblock_->asize_;

      // all slot in
      if (rank_ == 3) {
        double* const data = dfblock_->data_.get();
        for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {
          for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1) {
            for (int j2 = offset_[2]; j2 != offset_[2] + shell_[1]->nbasis(); ++j2, ++ppt) {
              data[j2+naux*(j1+nbin*j0)] = data[j2+naux*(j0+nbin*j1)] = *ppt;
            }
          }
        }
#if 0
      } else if (rank_ == 2) {
        double* const data = df_->data2_.get();
        for (int j0 = offset_[0]; j0 != offset_[0] + shell_[2]->nbasis(); ++j0) {
          for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++ppt) {
            data[j1+j0*naux] = data[j0+j1*naux] = *ppt;
          }
        }
#endif
      } else {
        assert(false);
      }
    };
};

}

#endif

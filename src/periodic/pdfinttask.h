//
// BAGEL - Parallel electron correlation program.
// Filename: pdfinttask.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#ifndef __SRC_PERIODIC_PDFINTTASK_H
#define __SRC_PERIODIC_PDFINTTASK_H

#include <src/df/dfblock.h>
#include <src/molecule/shell.h>
#include <src/integral/rys/eribatch.h>

namespace bagel {

class PDFIntTask_2index {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<int,2> offset_;
    DFDist* df_;

  public:
    // (i.|jL.) sum over L
    PDFIntTask_2index(std::array<std::shared_ptr<const Shell>,4>&& sh, std::array<int,2>&& offset, DFDist* df)
     : shell_(sh), offset_(offset), df_(df) { }

    void compute() {

      auto eribatch = std::make_shared<ERIBatch>(shell_, 2.0);
      const double* eridata = eribatch->data();

      const size_t naux = df_->naux();

      assert(offset_.size() == 3);
      double* const data = df_->data2_->data();
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[2]->nbasis(); ++j0)
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++eridata)
          data[j0+j1*naux] += *eridata;
    }
};


class PDFIntTask_3index {
  protected:
    const std::array<std::shared_ptr<const Shell>,4> shell_;
    const std::array<int,3> offset_;
    std::shared_ptr<DFBlock> dfblock_;

  public:
    // (r sL'|iL .) sum over L
    PDFIntTask_3index(std::array<std::shared_ptr<const Shell>,4>&& shells, std::array<int,3>&& offset, std::shared_ptr<DFBlock>& df)
     : shell_(shells), offset_(offset), dfblock_(df) { };

    void compute() {

      auto eribatch = std::make_shared<ERIBatch>(shell_, 2.0);

      assert(dfblock_->b1size() == dfblock_->b2size());
      const size_t nbin = dfblock_->b1size();
      const size_t naux = dfblock_->asize();
      const double* eridata = eribatch->data();

      double* const data = dfblock_->data();
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0)
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1)
          for (int j2 = offset_[2]; j2 != offset_[2] + shell_[1]->nbasis(); ++j2, ++eridata)
            data[j2 + naux * (j1 + nbin * j0)] += *eridata;
    }
};


}

#endif

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

class PDFIntTask {
  protected:
    const std::array<std::shared_ptr<const Shell>,4> shell_;
    const std::array<int,3> offset_; // at most 3 elements
    std::shared_ptr<DFBlock> dfblocks_;

    std::shared_ptr<ERIBatch> compute_batch(const std::array<std::shared_ptr<const Shell>,4>& input) const {
      auto eribatch = std::make_shared<ERIBatch>(input, 2.0);
      eribatch->compute();
      return eribatch;
    }

  public:
    PDFIntTask(std::array<std::shared_ptr<const Shell>,4>&& shells, std::array<int,3>&& offset, std::shared_ptr<DFBlock>& df)
     : shell_(shells), offset_(offset), dfblocks_(df) { };

    void compute() {
      std::shared_ptr<ERIBatch> eribatch = compute_batch(shell_);

      /* b1size = nbasis in cell 0 ** b2size = nbasis in cell g */
      assert(dfblocks_->b1size() == dfblocks_->b2size()); // problem here. Need to sum up the results.
      const size_t nbin = dfblocks_->b1size();
      const size_t naux = dfblocks_->asize();
      const double* eridata = eribatch->data();
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0)
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1, eridata += shell_[1]->nbasis())
          std::copy_n(eridata, shell_[1]->nbasis(), dfblocks_->data() + offset_[2] + naux * (j1 + nbin * j0));
    }

};


}

#endif

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
// The BAGEL package is free software; you can redistribute it and/or modify
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

#include <src/df/dfblock.h>
#include <src/scf/shell.h>
#include <src/rysint/rysint.h>

namespace bagel {

class DFIntTask_base {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<int,3> offset_; // at most 3 elements
    int rank_;
    std::vector<std::shared_ptr<DFBlock> > dfblocks_;

    virtual std::shared_ptr<Integral> compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) = 0;
    virtual int nblocks() const = 0;

  public:
    DFIntTask_base(std::array<std::shared_ptr<const Shell>,4>& a, std::vector<int>& b, std::vector<std::shared_ptr<DFBlock> >& df);

    void compute();

};

template <typename TBatch>
class DFIntTask : public DFIntTask_base {
  protected:
    std::shared_ptr<Integral> compute_batch(std::array<std::shared_ptr<const Shell>,4>& input) override {
      std::shared_ptr<TBatch> eribatch(new TBatch(input, 2.0));
      eribatch->compute();
      return eribatch;
    }

    int nblocks() const override { return TBatch::nblocks(); }

  public:
    DFIntTask(std::array<std::shared_ptr<const Shell>,4>&& a, std::vector<int>&& b, std::vector<std::shared_ptr<DFBlock> >& df) : DFIntTask_base(a,b,df) {} 

};


}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: mixederibatch.h
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


#ifndef __SRC_REL_MIXEDERIBATCH_H
#define __SRC_REL_MIXEDERIBATCH_H

#include <stddef.h>
#include <src/wfn/shell.h>
#include <src/wfn/geometry.h>
#include <memory>
#include <src/rysint/eribatch.h>
#include <src/util/matrix.h>
#include <src/parallel/resources.h>
#include <src/rysint/integral.h>

namespace bagel {

class MixedERIBatch : public Integral {
  protected:
    double* data_;

    // size info
    size_t size_block_;
    size_t size_alloc_;

    const std::array<std::shared_ptr<const Shell>,3> shells_;

    std::shared_ptr<StackMem> stack_;

  public:
    MixedERIBatch(std::array<std::shared_ptr<const Shell>,4> info, const double dummy);

    ~MixedERIBatch();

    void compute();

    double* data(const int i) override { assert(i < 3); return data_+i*size_block_; }

    void eri_compute(double* eri) const;

    size_t size_block() const { return size_block_; }
    size_t size_alloc() const { return size_alloc_; }

    // need to return 6 blocks
    constexpr static int nblocks() { return 3; }
};

}

#endif

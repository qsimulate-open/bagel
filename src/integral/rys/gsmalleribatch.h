//
// BAGEL - Parallel electron correlation program.
// Filename: gsmalleribatch.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_INTEGRAL_RYS_GSMALLERIBATCH_H
#define __SRC_INTEGRAL_RYS_GSMALLERIBATCH_H

#include <src/molecule/shell.h>
#include <src/math/xyzfile.h>
#include <src/parallel/resources.h>

namespace bagel {

class GSmallERIBatch {
  protected:
    double* data_;

    // size info
    size_t size_block_;
    size_t size_alloc_;

    // input shells
    const std::array<std::shared_ptr<const Shell>,3> shells_;

    // target atoms
    const std::array<int,3> atoms_;
    const int natoms_;

    std::shared_ptr<StackMem> stack_;

    double* data(const int i) { return data_+i*size_block_; }
    const double* data(const int i) const { return data_+i*size_block_; }

  public:
    GSmallERIBatch(std::array<std::shared_ptr<const Shell>,4> info, std::array<int,3>, const int);
    ~GSmallERIBatch();

    void compute();
    std::shared_ptr<GradFile> compute_gradient(std::array<std::shared_ptr<const Matrix>,6>&) const;

    size_t size_block() const { return size_block_; }
    constexpr static int nblocks() { return 9; }
};

}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: sphmmbatch.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __SRC_INTEGRAL_OS_SPHMMBATCH_H
#define __SRC_INTEGRAL_OS_SPHMMBATCH_H

#include <src/integral/os/mmbatch.h>
#include <src/util/parallel/resources.h>

namespace bagel {

class SphMMBatch {
  protected:
    std::array<std::shared_ptr<const Shell>,2> basisinfo_;
    std::shared_ptr<const Molecule> mol_;
    const int lmax_;

    double* data_;
    double* stack_save_;
    std::shared_ptr<StackMem> stack_;

    int nblocks() const { return (lmax_ + 1) * (lmax_ + 1) - 1; } //excluding monopoles
    size_t size_block_;
    size_t size_alloc_;

  public:
    SphMMBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, std::shared_ptr<const Molecule> m, const int lmax);
    ~SphMMBatch();

    void compute();

    int num_blocks() const { return nblocks(); }
    size_t size_block() const { return size_block_; }
    double* data(const int i) { assert(i < nblocks()); return data_+i*size_block_; }
    const double* data() const { return data_; }


};

class SphDipoleBatch : public SphMMBatch {
  public:
    SphDipoleBatch(const std::array<std::shared_ptr<const Shell>,2>& basis, std::shared_ptr<const Molecule> c) : SphMMBatch(basis, c, 1)
      { assert(nblocks() == Nblocks()); }
    constexpr static int Nblocks() { return 3; }
};

}

#endif

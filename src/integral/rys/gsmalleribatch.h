//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gsmalleribatch.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_INTEGRAL_RYS_GSMALLERIBATCH_H
#define __SRC_INTEGRAL_RYS_GSMALLERIBATCH_H

#include <src/molecule/shell.h>
#include <src/util/math/xyzfile.h>
#include <src/util/parallel/resources.h>
#include <src/util/math/btas_interface.h>

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
    std::shared_ptr<GradFile> compute_gradient(std::array<std::shared_ptr<const btas::Tensor3<double>>,6>&) const;

    size_t size_block() const { return size_block_; }
    constexpr static int nblocks() { return 9; }
};

}

#endif

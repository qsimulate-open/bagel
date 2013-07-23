//
// BAGEL - Parallel electron correlation program.
// Filename: smallnaibatch.h
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


#ifndef __SRC_INTEGRAL_RYS_SMALLNAIBATCH_H
#define __SRC_INTEGRAL_RYS_SMALLNAIBATCH_H

#include <stddef.h>
#include <src/molecule/shell.h>
#include <src/molecule/molecule.h>
#include <memory>
#include <src/integral/rys/naibatch.h>
#include <src/math/matrix.h>

// computes (sigma p)Vnuc(sigma p), and returns 4 blocks of data

namespace bagel {

class SmallNAIBatch {
  protected:
    std::array<std::shared_ptr<Matrix>,4> data_;

    const std::shared_ptr<const Molecule> mol_;
    const std::array<std::shared_ptr<const Shell>,2> shells_;

    const size_t size_block_;

    virtual std::shared_ptr<Matrix> nai_compute() const;

  public:
    SmallNAIBatch(std::array<std::shared_ptr<const Shell>,2> info, std::shared_ptr<const Molecule> mol);

    virtual void compute();

    std::shared_ptr<Matrix> operator[](const int i) { return data_[i]; }

    size_t size_block() const { return size_block_; }

};

}

#endif

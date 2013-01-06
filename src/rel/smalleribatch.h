//
// BAGEL - Parallel electron correlation program.
// Filename: smalleribatch.h
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


#ifndef __SRC_REL_SMALLERIBATCH_H
#define __SRC_REL_SMALLERIBATCH_H

#include <stddef.h>
#include <src/scf/shell.h>
#include <src/scf/geometry.h>
#include <memory>
#include <src/rysint/eribatch.h>
#include <src/rel/relshell.h>
#include <src/util/matrix.h>
#include <src/parallel/resources.h>

namespace bagel {

class SmallERIBatch {
  protected:
    double* data_;

    // size info
    size_t size_block_;
    size_t size_alloc_;

    const std::array<std::shared_ptr<const RelShell>,4> shells_;

    std::shared_ptr<StackMem> stack_;

  public:
    SmallERIBatch(std::array<std::shared_ptr<const RelShell>,4> info);

    ~SmallERIBatch();

    void compute();

    void eri_compute(double* eri) const;

    size_t size_block() const { return size_block_; }
    size_t size_alloc() const { return size_alloc_; }
};

}

#endif

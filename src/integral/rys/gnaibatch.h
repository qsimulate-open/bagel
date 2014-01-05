//
// BAGEL - Parallel electron correlation program.
// Filename: gnaibatch.h
// Copyright (C) 2009 Toru Shiozaki
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

#ifndef __SRC_GRAD_GNAIBATCH_H
#define __SRC_GRAD_GNAIBATCH_H

#include <memory>
#include <tuple>
#include <src/math/xyzfile.h>
#include <src/integral/rys/coulombbatch_base.h>

namespace bagel {

class GNAIBatch : public CoulombBatch_base {

  protected:
    void set_exponents();
    std::unique_ptr<double[]> exponents_;

    std::tuple<int,int> iatom_;

  public:

    GNAIBatch(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, const std::tuple<int,int> i,
              std::shared_ptr<StackMem> stack = std::shared_ptr<StackMem>());

    /// compute a batch of integrals
    void compute();

    std::shared_ptr<GradFile> compute_gradient(std::shared_ptr<const Matrix> cden, const int iatom0, const int iatom1, const int natom) const;

    int nblocks() const { return mol_->natom()*3; }

    double* data(const int i) override { return data_ + i*size_block_; }
};

}

#endif


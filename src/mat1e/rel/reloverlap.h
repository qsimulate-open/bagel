//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reloverlap.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_REL_RELOVERLAP_H
#define __SRC_REL_RELOVERLAP_H

#include <src/mat1e/overlap.h>
#include <src/mat1e/kinetic.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

class RelOverlap : public ZMatrix {
  protected:
    const std::shared_ptr<const Molecule> mol_;
    const std::shared_ptr<const Matrix> kinetic_;
    const std::shared_ptr<const Overlap> overlap_;

    void compute_();

  public:
    RelOverlap(const std::shared_ptr<const Molecule> mol) : ZMatrix(mol->nbasis()*4, mol->nbasis()*4),
               mol_(mol), kinetic_(std::make_shared<Kinetic>(mol)), overlap_(std::make_shared<Overlap>(mol)) {
      compute_();
    }

    std::shared_ptr<ZMatrix> tildex(const double thresh) const override;
    void inverse() override;
};

}

#endif

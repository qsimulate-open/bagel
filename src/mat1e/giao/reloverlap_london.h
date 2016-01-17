//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reloverlap_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_LONDON_RELOVERLAP_LONDON_H
#define __SRC_LONDON_RELOVERLAP_LONDON_H

#include <src/mat1e/giao/zoverlap.h>
#include <src/mat1e/giao/zkinetic.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

class RelOverlap_London : public ZMatrix {
  protected:
    const std::shared_ptr<const Molecule> mol_;
    const std::shared_ptr<const ZMatrix> kinetic_;
    const std::shared_ptr<const ZOverlap> overlap_;

    void compute_();

  public:
    RelOverlap_London(const std::shared_ptr<const Molecule> mol) : ZMatrix(mol->nbasis()*4, mol->nbasis()*4),
                      mol_(mol), kinetic_(std::make_shared<ZKinetic>(mol)), overlap_(std::make_shared<ZOverlap>(mol)) {
      compute_();
    }

    std::shared_ptr<ZMatrix> tildex(const double thresh) const;
    void inverse() override;
};

}

#endif

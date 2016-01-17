//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relhcore_london.h
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


#ifndef __SRC_LONDON_RELHCORE_LONDON_H
#define __SRC_LONDON_RELHCORE_LONDON_H

#include <src/mat1e/giao/zhcore.h>
#include <src/mat1e/giao/zkinetic.h>
#include <src/mat1e/giao/zoverlap.h>
#include <src/mat1e/giao/small1e_london.h>
#include <src/integral/comprys/complexnaibatch.h>
#include <src/integral/compos/complexoverlapbatch.h>

namespace bagel {

class RelHcore_London : public ZMatrix {
  protected:
    const std::shared_ptr<const Molecule> geom_;
    const std::shared_ptr<const ZMatrix> kinetic_;
    const std::shared_ptr<const ZMatrix> hcore_;
    const std::shared_ptr<const ZMatrix> nai_;
    const std::shared_ptr<const ZMatrix> overlap_;
    std::shared_ptr<Small1e_London<ComplexNAIBatch>> smallnai_;

    void compute_();


  public:
    RelHcore_London(const std::shared_ptr<const Molecule> geom) : ZMatrix(geom->nbasis()*4, geom->nbasis()*4), geom_(geom),
            kinetic_(std::make_shared<ZKinetic>(geom_)),
            hcore_(std::make_shared<ZHcore>(geom_)),
            nai_(std::make_shared<ZMatrix>(*hcore_ - *kinetic_)),
            overlap_(std::make_shared<ZOverlap>(geom_)) {
      smallnai_ = std::make_shared<Small1e_London<ComplexNAIBatch>>(geom_);
      if (geom_->has_finite_nucleus()) {
        smallnai_->ax_plus_y(1.0, *std::make_shared<Small1e_London<ComplexERIBatch>>(geom_));
      }
      compute_();
    }

};

}

#endif

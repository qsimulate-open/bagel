//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relhcore.h
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


#ifndef __SRC_REL_RELHCORE_H
#define __SRC_REL_RELHCORE_H

#include <src/mat1e/hcore.h>
#include <src/mat1e/kinetic.h>
#include <src/mat1e/nai.h>
#include <src/mat1e/rel/small1e.h>
#include <src/integral/rys/naibatch.h>

namespace bagel {

class RelHcore : public ZMatrix {
  protected:
    const std::shared_ptr<const Molecule> geom_;
    const std::shared_ptr<const Matrix> kinetic_;
    const std::shared_ptr<const Matrix> nai_;
    std::shared_ptr<Small1e<NAIBatch>> smallnai_;

    void compute_();

  public:
    RelHcore(std::shared_ptr<const Molecule> geom)
     : ZMatrix(geom->nbasis()*4, geom->nbasis()*4), geom_(geom),
       kinetic_(std::make_shared<Kinetic>(geom_)), nai_(std::make_shared<NAI>(geom_)) {
      smallnai_ = std::make_shared<Small1e<NAIBatch>>(geom_);
      if (geom_->has_finite_nucleus())
        smallnai_->ax_plus_y(1.0, *std::make_shared<Small1e<ERIBatch>>(geom_));
      compute_();
    }

};

}

#endif

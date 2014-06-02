//
// BAGEL - Parallel electron correlation program.
// Filename: relhcore_london.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_LONDON_RELHCORE_LONDON_H
#define __SRC_LONDON_RELHCORE_LONDON_H

#include <src/molecule/zhcore.h>
#include <src/molecule/zkinetic.h>
#include <src/molecule/zoverlap.h>
#include <src/integral/comprys/complexnaibatch.h>
#include <src/integral/compos/complexoverlapbatch.h>
#include <src/london/small1e_london.h>

namespace bagel {

class RelHcore_London : public ZMatrix {
  protected:
    const std::shared_ptr<const Molecule> geom_;
    const std::shared_ptr<const ZMatrix> kinetic_;
    const std::shared_ptr<const ZMatrix> hcore_;
    const std::shared_ptr<const ZMatrix> nai_;
    const std::shared_ptr<const ZMatrix> overlap_;
    std::shared_ptr<Small1e_London<ComplexNAIBatch>> smallnai_;
    std::shared_ptr<Small1e_London<ComplexOverlapBatch>> smalloverlap_;

    void compute_();


  public:
    RelHcore_London(const std::shared_ptr<const Molecule> geom) : ZMatrix(geom->nbasis()*4, geom->nbasis()*4), geom_(geom),
            kinetic_(std::make_shared<ZKinetic>(geom_)),
            hcore_(std::make_shared<ZHcore>(geom_)),
            nai_(std::make_shared<ZMatrix>(*hcore_ - *kinetic_)),
            overlap_(std::make_shared<ZOverlap>(geom_)) {
      smallnai_ = std::make_shared<Small1e_London<ComplexNAIBatch>>(geom_, false);
      smalloverlap_ = std::make_shared<Small1e_London<ComplexOverlapBatch>>(geom_, true);
      if (geom_->has_finite_nucleus()) {
        smallnai_->ax_plus_y(1.0, *std::make_shared<Small1e_London<ComplexERIBatch>>(geom_, false));
      }
      compute_();
    }

};

}

#endif

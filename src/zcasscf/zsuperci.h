//
// BAGEL - Parallel electron correlation program.
// Filename: zsuperci.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson Bates <jefferson.bates@northwestern.edu>
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

#ifndef __SRC_ZCASSCF_ZSUPERCI_H
#define __SRC_ZCASSCF_ZSUPERCI_H

#include <src/zcasscf/zcasscf.h>

namespace bagel {

class ZSuperCI : public ZCASSCF {
  protected:
    void common_init() {
      std::cout << "   * Using the Super CI algorithm as noted in Roos (1980) IJQC" << std::endl << std::endl;
    }

    // internal function
    // gradients

    // diagonal denominator

  public:
    ZSuperCI(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> ref = nullptr) 
       : ZCASSCF(idat, geom, ref) { common_init(); }

    void compute() override;
    void one_body_operators(std::shared_ptr<ZMatrix>& f, std::shared_ptr<ZMatrix>& fact, std::shared_ptr<ZMatrix>& factp, std::shared_ptr<ZMatrix>& gaa,
                            std::shared_ptr<ZRotFile>& denom);

};

}

#endif

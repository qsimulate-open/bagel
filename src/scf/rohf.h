//
// BAGEL - Parallel electron correlation program.
// Filename: rohf.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#ifndef __BAGEL_SRC_SCF_ROHF_H
#define __BAGEL_SRC_SCF_ROHF_H

#include <src/scf/uhf.h>

// implements UHF as in Tsuchimochi and Scuseria, J. Chem. Phys. 133, 141102 (2010)
namespace bagel {

class ROHF : public UHF {
  protected:
    void symmetrize_cv(std::shared_ptr<Matrix>, std::shared_ptr<Matrix>);

  public:
    ROHF(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom,
         const std::shared_ptr<const Reference> re = std::shared_ptr<const Reference>()) : UHF(idata, geom, re) {
      // TODO projection does not work
      coeff_ = std::shared_ptr<Coeff>();
      coeffB_ = std::shared_ptr<Coeff>();
    }

    void compute() override;
};

}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbopt/ras/bfgs.h
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Inkoo Kim: <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __BAGEL_ASD_ORBOPT_BFGS_H
#define __BAGEL_ASD_ORBOPT_BFGS_H

#include <src/asd/orbopt/ras/rasscf.h>

namespace bagel {

class ASDRASBFGS : public ASDRASSCF {

  protected:
    void common_init() {
      std::cout << "    * Using the Quasi 2nd-order algorithm as noted in Chaban et al. TCA (1997)" << std::endl << std::endl;
    }

    // compute orbital gradients
    void grad_vc(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<ASDRASRotFile> sigma) const;
    void grad_va(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> qxr,   std::shared_ptr<Matrix> rdm1, std::shared_ptr<ASDRASRotFile> sigma) const;
    void grad_ca(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<Matrix> rdm1, std::shared_ptr<ASDRASRotFile> sigma) const;

    void grad_aa(std::shared_ptr<const Matrix> mcfock, std::shared_ptr<ASDRASRotFile> sigma) const;

    // compute diagonal denominators
    std::shared_ptr<const ASDRASRotFile> compute_denom(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<const Matrix> rdm1, std::shared_ptr<const Matrix> mcfock) const;

  public:
    ASDRASBFGS(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
      : ASDRASSCF(idat, geom, ref) { 
      std::cout << "ASDRASBFGS(active-active) constructor" << std::endl; 
      common_init(); 
    }

    void compute() override;

};

}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: asd/orbital/rasbfgs.h
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

#ifndef __BAGEL_ASD_RAS_BFGS_H
#define __BAGEL_ASD_RAS_BFGS_H

#include <src/asd/orbital/rasscf.h>

namespace bagel {

class ASD_RAS_BFGS : public ASD_RASSCF {

  protected:
    void common_init() {
      std::cout << "    * Using the Quasi 2nd-order algorithm as noted in Chaban et al. TCA (1997)" << std::endl << std::endl;
    }

    // compute orbital gradients
    void grad_vc(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<ASD_RAS_RotFile> sigma) const;
    void grad_va(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> qxr,   std::shared_ptr<Matrix> rdm1, std::shared_ptr<ASD_RAS_RotFile> sigma) const;
    void grad_ca(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<Matrix> rdm1, std::shared_ptr<ASD_RAS_RotFile> sigma) const;

    void grad_aa(std::shared_ptr<const Matrix> mcfock, std::shared_ptr<ASD_RAS_RotFile> sigma) const;

    //large
    void grad_vc_large(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<RotFile> sigma) const;
    void grad_va_large(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> qxr,   std::shared_ptr<Matrix> rdm1, std::shared_ptr<RotFile> sigma) const;
    void grad_ca_large(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<Matrix> rdm1, std::shared_ptr<RotFile> sigma) const;

    //small
    void grad_aa_small(std::shared_ptr<const Matrix> mcfock, std::shared_ptr<ASD_RAS_ActiveRotFile> sigma) const;

    // compute diagonal denominators
    std::shared_ptr<const ASD_RAS_RotFile> compute_denom(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<const Matrix> rdm1, std::shared_ptr<const Matrix> mcfock) const;

    std::shared_ptr<const RotFile> compute_denom_large(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<const Matrix> rdm1) const;
    std::shared_ptr<const ASD_RAS_ActiveRotFile> compute_denom_small(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> rdm1, std::shared_ptr<const Matrix> mcfock) const;
  public:
    ASD_RAS_BFGS(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
      : ASD_RASSCF(idat, geom, ref) { 
      std::cout << "ASD_RAS_BFGS(active-active) constructor" << std::endl; 
      common_init(); 
    }

    void compute() override;

};

}

#endif

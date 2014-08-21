//
// BAGEL - Parallel electron correlation program.
// Filename: zcasbfgs.h
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

#ifndef __SRC_ZCASSCF_ZCASBFGS_H
#define __SRC_ZCASSCF_ZCASBFGS_H

#include <src/zcasscf/zcasscf.h>

namespace bagel {

class ZCASBFGS : public ZCASSCF {
  protected:
    void common_init() {
      std::cout << "   * Using the Quasi 2nd-order algorithm as noted in Chaban et al. TCA (1997)" << std::endl << std::endl;
    }

    // internal function
    // gradient functions for computing the virtual-closed, virtual-active, and closed-active contributions
    void grad_vc(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<ZRotFile> sigma) const;
    void grad_va(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> qxr,   std::shared_ptr<const ZMatrix> rdm1, std::shared_ptr<ZRotFile> sigma) const;
    void grad_ca(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<const ZMatrix> qxr,
                 std::shared_ptr<const ZMatrix> rdm1, std::shared_ptr<ZRotFile> sigma) const;

    // diagonal Hessian
    std::shared_ptr<ZRotFile> compute_denom(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock,
                 std::shared_ptr<const ZMatrix> qxr, std::shared_ptr<const ZMatrix> rdm1) const;

  public:
    ZCASBFGS(const std::shared_ptr<const PTree> idat, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> ref = nullptr)
       : ZCASSCF(idat, geom, ref) { common_init(); }

    void compute() override;
   // function to optimize only the electronic-electronic or electronic-positronic type rotations
   std::tuple<std::shared_ptr<ZRotFile>, std::vector<double>, std::shared_ptr<ZRotFile>, std::shared_ptr<ZRotFile>, bool> optimize_subspace_rotations(std::vector<double> energy, std::shared_ptr<const ZRotFile> grad, std::shared_ptr<const ZRotFile> rot, std::shared_ptr<SRBFGS<ZRotFile>> srbfgs, std::shared_ptr<ZMatrix> cold, bool optimize_electrons = true);
   // function to copy electronic rotations from a rotation file TODO: make lambda
   std::shared_ptr<ZRotFile> copy_electronic_rotations(std::shared_ptr<const ZRotFile> rot) const;
   // function to copy positronic rotations from a rotation file TODO: make lambda
   std::shared_ptr<ZRotFile> copy_positronic_rotations(std::shared_ptr<const ZRotFile> rot) const;
   // returns "optimal" level shift
   std::complex<double> find_level_shift(std::shared_ptr<const ZRotFile> rotmat) const;

};

}

#endif

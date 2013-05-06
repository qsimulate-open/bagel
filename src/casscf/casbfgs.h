//
// BAGEL - Parallel electron correlation program.
// Filename: casbfgs.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __BAGEL_CASSCF_CASBFGS_H
#define __BAGEL_CASSCF_CASBFGS_H

#include <src/casscf/casscf.h>

namespace bagel {

class CASBFGS : public CASSCF {

  protected:
    void common_init() {
      std::cout << "    * Using the Quasi 2nd-order algorithm as noted in Chaban et al. TCA (1997)" << std::endl;
    }

    // compute orbital gradients
    void grad_vc(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<RotFile> sigma) const;
    void grad_va(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> qxr,   std::shared_ptr<RotFile> sigma) const;
    void grad_ca(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr, std::shared_ptr<RotFile> sigma) const;

    // compute diagonal denominators
    std::shared_ptr<const RotFile> compute_denom(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr) const;

  public:
    CASBFGS(const boost::property_tree::ptree& idat, const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const Reference> ref)
      : CASSCF(idat, geom, ref) { common_init(); }

    void compute();

};

}

#endif

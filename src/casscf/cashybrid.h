//
// BAGEL - Parallel electron correlation program.
// Filename: casbfgs.h
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


#ifndef __BAGEL_CASSCF_CASHYBRID_H
#define __BAGEL_CASSCF_CASHYBRID_H

#include <iostream>
#include <iomanip>
#include <src/wfn/method.h>

namespace bagel {

class CASHYBRID : public Method {

  protected:
    int maxiter_sci_;
    int maxiter_bfgs_;
    double thresh_switch_;

    std::shared_ptr<const Reference> refout_;

    void common_init() {
      std::cout << "    * Using a hybrid approach to CASSCF *    " << std::endl;
      maxiter_sci_   = idata_->get<int>("maxiter_sci", 10);
      maxiter_bfgs_  = idata_->get<int>("maxiter_bfgs", 10);
      thresh_switch_ = idata_->get<double>("thresh_switch", 1.0e-4);
      std::cout << "    * SuperCI used for a maximum of " << maxiter_sci_ << " iterations *    " << std::endl;
      std::cout << std::setprecision(2) << std::scientific
                << "    * Switching to qusi-Newton methods when the rms gradient norm is smaller than " << thresh_switch_ << " *    " << std::endl;
      std::cout << "    * Step-restricted BFGS will then be used for a maximum of " << maxiter_bfgs_ << " iterations *    " << std::endl << std::endl;
    }


  public:
    CASHYBRID(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
      : Method(idat, geom, ref) { common_init(); }

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const;
};

}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: werner.h
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


#ifndef __BAGEL_CASSCF_WERNER_H
#define __BAGEL_CASSCF_WERNER_H

#include <src/casscf/casscf.h>
#include <src/casscf/rotfile.h>
#include <src/casscf/jvec.h>

namespace bagel {

class WernerKnowles : public CASSCF {
  protected:
    virtual void common_init() {
      std::cout << "    * Using the two-step Werner-Knowles algorithm (see JCP 1985)" << std::endl << std::endl;
    };

    std::shared_ptr<Matrix> compute_bvec(const std::shared_ptr<const Jvec>, std::shared_ptr<Matrix>, std::shared_ptr<const Coeff>);
    std::shared_ptr<Matrix> compute_bvec(const std::shared_ptr<const Jvec>, std::shared_ptr<Matrix>,
                                           std::shared_ptr<Matrix>, const std::shared_ptr<const Coeff>);
    std::shared_ptr<Matrix> compute_bvec2(std::shared_ptr<Jvec>, std::shared_ptr<Matrix>,
                                           std::shared_ptr<Matrix>, const std::shared_ptr<const Coeff>);
    std::shared_ptr<const Matrix> compute_denom(const std::shared_ptr<const Matrix>);
    std::shared_ptr<Matrix> compute_sigma_R(const std::shared_ptr<const Jvec>, const std::shared_ptr<const Matrix>,
                                              const std::shared_ptr<const Matrix>, const std::shared_ptr<const Matrix>);

    double thresh_mmicro_;
    int max_mmicro_iter_;

  public:
    WernerKnowles(const boost::property_tree::ptree& idat, const std::shared_ptr<const Geometry> geom)
      : CASSCF(idat, geom) {common_init();
      // get thresh (for micro iteration) from the input
      thresh_mmicro_ = idat.get<double>("thresh_mmicro", thresh_micro_);
      max_mmicro_iter_ = idat.get<int>("maxiter_mmicro", 3);
    };
    ~WernerKnowles() {};

    virtual void compute();


};

}

#endif


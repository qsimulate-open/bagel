//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcassecond.h
// Copyright (C) 2016 Toru Shiozaki
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

#ifndef __SRC_ZCASSCF_ZCASSecond_H
#define __SRC_ZCASSCF_ZCASSecond_H

#include <src/multi/zcasscf/zcasscf.h>

namespace bagel {

class ZCASSecond : public ZCASSCF {
  protected:
    // convergence threshold for micro iteration relative to stepsize
    double thresh_microstep_;

    // compute orbital gradient
    std::shared_ptr<ZRotFile> compute_gradient(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<const ZMatrix> qxr) const;

    // diagonal Hessian
    std::shared_ptr<ZRotFile> compute_denom(std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock,
                                            std::shared_ptr<const ZMatrix> qxr, std::shared_ptr<const ZMatrix> rdm1) const;
    // compute H*t (Hessian times trial vector)
    std::shared_ptr<ZRotFile> compute_hess_trial(std::shared_ptr<const ZRotFile> trot,
              std::list<std::shared_ptr<const RelDFHalf>> halfc, std::list<std::shared_ptr<const RelDFHalf>> halfg, std::list<std::shared_ptr<const RelDFHalf>> halfb,
              std::list<std::shared_ptr<const RelDFHalf>> halfac, std::list<std::shared_ptr<const RelDFHalf>> halfag, std::list<std::shared_ptr<const RelDFHalf>> halfab,
              std::shared_ptr<const ZMatrix> cfock, std::shared_ptr<const ZMatrix> afock, std::shared_ptr<const ZMatrix> qxr) const;
    // apply denominator in microiterations
    std::shared_ptr<ZRotFile> apply_denom(std::shared_ptr<const ZRotFile> grad, std::shared_ptr<const ZRotFile> denom, const double shift, const double scale) const;

  public:
    ZCASSecond(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
       : ZCASSCF(idat, geom, ref) {
      std::cout << "   * Using the second-order algorithm" << std::endl << std::endl;
      // overwriting thresh_micro
      thresh_micro_ = idata_->get<double>("thresh_micro", thresh_*0.5);
      thresh_microstep_ = idata_->get<double>("thresh_microstep", 1.0e-4);
    }

    void compute() override;

};

}

#endif

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cassecond.h
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


#ifndef __BAGEL_CASSCF_SECOND_H
#define __BAGEL_CASSCF_SECOND_H

#include <src/multi/casscf/casscf.h>
#include <src/wfn/rdm.h>

namespace bagel {

// implements the second-order algorithm with augmented Hessian (with the help of Takeshi Yanai)

class CASSecond : public CASSCF {
  protected:
    // convergence threshold for micro iteration relative to stepsize
    double thresh_microstep_;

    // compute orbital gradient
    std::shared_ptr<RotFile> compute_gradient(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr) const;
    // compute exact diagonal Hessian
    std::shared_ptr<RotFile> compute_denom(std::shared_ptr<const DFHalfDist> half, std::shared_ptr<const DFHalfDist> half_1j, std::shared_ptr<const DFHalfDist> halfa,
                                           std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock) const;
    // compute H*t (Hessian times trial vector)
    std::shared_ptr<RotFile> compute_hess_trial(std::shared_ptr<const RotFile> trot, std::shared_ptr<const DFHalfDist> half, std::shared_ptr<const DFHalfDist> halfa,
                                                std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr) const;
    // apply denominator in microiterations
    std::shared_ptr<RotFile> apply_denom(std::shared_ptr<const RotFile> grad, std::shared_ptr<const RotFile> denom, const double shift, const double scale) const;

  public:
    CASSecond(std::shared_ptr<const PTree> idat, std::shared_ptr<const Geometry> geom, std::shared_ptr<const Reference> ref = nullptr)
      : CASSCF(idat, geom, ref) {
      std::cout << "    * Using the second-order algorithm" << std::endl << std::endl;
      // overwriting thresh_micro
      thresh_micro_ = idata_->get<double>("thresh_micro", thresh_*0.5);
      thresh_microstep_ = idata_->get<double>("thresh_microstep", 1.0e-4);
    }

    void compute() override;

    void trans_natorb();
};

}

#endif

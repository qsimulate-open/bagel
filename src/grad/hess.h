//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hess.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_GRAD_HESS_H
#define __SRC_GRAD_HESS_H

#include <src/wfn/reference.h>
#include <src/util/muffle.h>
#include <src/util/atommap.h>

namespace bagel {

class Hess {
  protected:
    const std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const Coeff> coeff_;

    bool numhess_;
    bool numforce_;

    int nproc_;

    std::shared_ptr<Matrix> hess_;
    std::shared_ptr<Matrix> mw_hess_;
    std::shared_ptr<Matrix> proj_hess_;
    std::shared_ptr<Matrix> eigvec_cart_;
    std::shared_ptr<Matrix> cartesian_;
    std::vector<double> ir_;
    std::vector<double> freq_;

    double dx_;
    double energy_;

    // mask some of the output
    mutable std::shared_ptr<Muffle> muffle_;

    void compute_finite_diff_();
    void project_zero_freq_();
    void print_ir_() const;

  public:
    Hess(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    std::shared_ptr<const Matrix> hess() const { return hess_; }
    std::shared_ptr<const Matrix> mw_hess() const { return mw_hess_; }
    std::shared_ptr<const Matrix> proj_hess() const { return proj_hess_; }

    void compute();
    virtual std::shared_ptr<const Reference> conv_to_ref() const;

};


}

#endif

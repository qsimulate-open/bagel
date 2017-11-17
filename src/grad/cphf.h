//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cphf.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_GRAD_CPHF_H
#define __SRC_GRAD_CPHF_H

#include <src/wfn/reference.h>

namespace bagel {

class CPHF {
  protected:
    std::shared_ptr<const Matrix> grad_;
    VectorB eig_;
    std::shared_ptr<const DFHalfDist> halfjj_;
    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const Geometry> geom_;

  public:
    CPHF(const std::shared_ptr<const Matrix> grad, const VectorB& eig,
         const std::shared_ptr<const DFHalfDist> half, const std::shared_ptr<const Reference> g);

    std::shared_ptr<Matrix> solve(const double thresh, const int maxiter = 100);

};

}

#endif


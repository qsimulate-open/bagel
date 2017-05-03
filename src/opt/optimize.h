//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: optimize.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_OPT_OPTIMIZE_H
#define __SRC_OPT_OPTIMIZE_H

#include <src/wfn/reference.h>

namespace bagel {

class Optimize {
  protected:
    const std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;

    int maxiter_;
    double thresh_;

  public:
    Optimize(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute();
    std::shared_ptr<const Geometry> geometry() const { return geom_; }
    std::shared_ptr<const Reference> conv_to_ref() const { return ref_; }
};


}


#endif

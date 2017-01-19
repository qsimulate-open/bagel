//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: force.h
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

#ifndef __SRC_GRAD_FORCE_H
#define __SRC_GRAD_FORCE_H

#include <src/wfn/reference.h>

namespace bagel {

class Force {
  protected:
    const std::shared_ptr<const PTree> idata_;
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const Reference> ref_;
    bool numerical_;		// numerical or analytical gradient?

  public:
    Force(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    std::shared_ptr<GradFile> compute();

};


}

#endif

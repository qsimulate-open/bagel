//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dmp2grad.h
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


#ifndef __SRC_REL_DMP2GRAD_H
#define __SRC_REL_DMP2GRAD_H

#include <src/pt2/dmp2/dmp2.h>
#include <src/wfn/reference.h>

namespace bagel {

class DMP2Grad : public DMP2 {
  protected:

  public:
    DMP2Grad(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }
};

}

#endif

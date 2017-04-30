//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mp2grad.h
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


#ifndef __SRC_MP2_MP2GRAD_H
#define __SRC_MP2_MP2GRAD_H

#include <src/pt2/mp2/mp2.h>
#include <src/wfn/reference.h>

namespace bagel {

class MP2Grad : public MP2 {
  protected:
  std::vector<double> dipole_;

  public:
    MP2Grad(std::shared_ptr<const PTree>, std::shared_ptr<const Geometry>, std::shared_ptr<const Reference>);

    void compute() override;

    std::shared_ptr<const Reference> conv_to_ref() const override { return ref_; }

    const std::vector<double>& dipole() const { return dipole_; }
    double dipole(const int i) const { return dipole_[i]; }

};

}

#endif

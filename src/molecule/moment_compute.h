//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moment_compute.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#ifndef __SRC_MOLECULE_MOMENT_COMPUTE_H
#define __SRC_MOLECULE_MOMENT_COMPUTE_H

#include <array>
#include <src/molecule/shell.h>
#include <src/util/math/matrix.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

struct MomentCompute {
  private:
    static std::array<std::shared_ptr<const Matrix>,6> mblock(const Shell& shell, const double exponent);
    static std::array<std::shared_ptr<const ZMatrix>,9> mblock(const Shell& shell, const double exponent, const std::array<double,3> magnetic_field, const bool london);
  public:
    // Gaussian orbitals
    static std::array<std::shared_ptr<const Matrix>,3> call(const Shell& shell);
    // London orbitals (and common origin)
    static std::array<std::shared_ptr<const ZMatrix>,3> call(const Shell& shell, const std::array<double,3> magnetic_field, const bool london);
};

}

#endif

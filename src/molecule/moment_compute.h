//
// BAGEL - Parallel electron correlation program.
// Filename: moment_compute.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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


#ifndef __SRC_MOLECULE_MOMENT_COMPUTE_H
#define __SRC_MOLECULE_MOMENT_COMPUTE_H

#include <array>
#include <src/math/matrix.h>
#include <src/math/zmatrix.h>

namespace bagel {

class Shell;

// Gaussian orbitals
std::array<std::shared_ptr<const Matrix>,3> moment_compute(const Shell& shell);
std::array<std::shared_ptr<const Matrix>,6> mblock(const Shell& shell, const double exponent);

// London orbitals (and common origin)
std::array<std::shared_ptr<const ZMatrix>,3> moment_compute(const Shell& shell, const std::array<double,3> magnetic_field, const bool london);
std::array<std::shared_ptr<const ZMatrix>,9> mblock(const Shell& shell, const double exponent, const std::array<double,3> magnetic_field, const bool london);

}

#endif

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhgrad.h
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#ifndef __SRC_GRAD_DKHGRAD_H
#define __SRC_GRAD_DKHGRAD_H

#include <src/molecule/molecule.h>
#include <src/util/math/diagvec.h>

namespace bagel {

class DKHgrad {
  private:
    // Number of uncontracted basis functions
    int nbasis_;

    // Transformation matrix to momentum space
    Matrix wtrans_;
    // Reverse transformation from momentum space
    Matrix wtrans_rev_;
    // Projection operator from uncontracted to contracted AO basis
    Matrix ptrans_;

    // Integrals in momentum space
    DiagVec kinetic_;
    Matrix nai_;
    Matrix smallnai_;

    // Lagrange multiplier for diagonalization of kinetic operator
    Matrix zmult_;
    // Energy derivative with respect to momentum orbital rotation
    Matrix ederiv_;

  public:
    DKHgrad() { }
    DKHgrad(std::shared_ptr<const Molecule>);

    // Effective density matrix for kinetic gradient
    std::shared_ptr<const Matrix> compute_tden(std::shared_ptr<const Matrix>);
    // Effective density matrices for NAI/SmallNAI gradients
    std::array<std::shared_ptr<const Matrix>, 2> compute_vden(std::shared_ptr<const Matrix>);
    // Effective density matrix for overlap gradient
    std::shared_ptr<const Matrix> compute_sden(std::shared_ptr<const Matrix>, std::shared_ptr<const Matrix>);
};

}

#endif

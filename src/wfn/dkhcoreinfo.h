//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcoreinfo.h
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


#ifndef __SRC_WFN_DKHCOREINFO_H
#define __SRC_WFN_DKHCOREINFO_H

#include <src/molecule/molecule.h>

namespace bagel {

class DKHcoreInfo {
  private:
    int nbasis_;

    Matrix wtrans_;
    Matrix ptrans_;

    VectorB kinetic_;
    Matrix nai_;
    Matrix smallnai_;

    Matrix zmult_;
    Matrix ederiv_;

  public:
    DKHcoreInfo() { }
    DKHcoreInfo(std::shared_ptr<const Molecule>);

    std::shared_ptr<const Matrix> compute_tden(std::shared_ptr<const Matrix>);
    std::array<std::shared_ptr<const Matrix>, 2> compute_vden(std::shared_ptr<const Matrix>);
    std::shared_ptr<const Matrix> compute_sden(std::shared_ptr<const Matrix>);
};

}

#endif

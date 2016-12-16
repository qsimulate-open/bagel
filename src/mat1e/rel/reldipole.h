//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reldipole.h
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

#ifndef __SRC_REL_RELDIPOLE_H
#define __SRC_REL_RELDIPOLE_H

#include <string>
#include <src/molecule/molecule.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

class RelDipole {
  protected:
    std::shared_ptr<const Molecule> geom_;
    std::shared_ptr<const ZMatrix> density_;
    std::string jobname_;

  public:
    RelDipole(std::shared_ptr<const Molecule> g, std::shared_ptr<const ZMatrix> z = nullptr, const std::string jobname = "") : geom_(g), density_(z), jobname_(jobname) { }

    std::array<double,3> compute() const;
    std::array<std::shared_ptr<const ZMatrix>,3> compute_matrices() const;
};

}

#endif

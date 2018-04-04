//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dfpcmo.h
// Copyright (C) 2018 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_UTIL_IO_DFPCMO_H
#define __SRC_UTIL_IO_DFPCMO_H

#include <src/util/math/zmatrix.h>

namespace bagel {

class DFPCMO {
  protected:
    std::shared_ptr<const ZMatrix> coeff_;
    std::shared_ptr<const VectorB> eig_;
    double energy_;
    int nesol_;
    int npsol_;
    int nbasis_;

  public:
    DFPCMO(std::shared_ptr<const ZMatrix> c, std::shared_ptr<const VectorB> e, double en, int ne, int np, int nb);
    void print() const;
};

}

#endif

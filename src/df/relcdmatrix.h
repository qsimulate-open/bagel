//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relcdmatrix.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew Kelley <matthewkelley2017@northwestern.edu>
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


#ifndef __SRC_DF_CDMATRIX_H
#define __SRC_DF_CDMATRIX_H

#include <src/df/reldfhalf.h>
#include <src/util/math/vectorb.h>

namespace bagel {

class RelDFHalf;

class RelCDMatrix : public ZVectorB {
  protected:
    const int alpha_comp_;

  public:
    RelCDMatrix(std::shared_ptr<const RelDFHalf> dfhc, std::shared_ptr<const SpinorInfo> abc, std::array<std::shared_ptr<const Matrix>, 4> trcoeff,
                std::array<std::shared_ptr<const Matrix>, 4> ticoeff, std::shared_ptr<const Matrix> dat2, const int number_of_j);
    RelCDMatrix(const ZVectorB& o, const int acomp) : ZVectorB(o), alpha_comp_(acomp) { }

    int alpha_comp() const { return alpha_comp_; }

};

}

#endif

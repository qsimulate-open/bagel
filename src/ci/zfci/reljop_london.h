//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: reljop_london.h
// Copyright (C) 2017 Toru Shiozaki
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

#ifndef __BAGEL_CI_ZFCI_RELJOP_LONDON_H
#define __BAGEL_CI_ZFCI_RELJOP_LONDON_H

#include <src/ci/zfci/reljop.h>
#include <src/mat1e/giao/relhcore_london.h>

namespace bagel {

class RelJop_London : public RelJop {
  protected:
    std::shared_ptr<ZMatrix> compute_hcore() const override { return std::make_shared<RelHcore_London>(geom_); }

  public:
    RelJop_London(const std::shared_ptr<const Geometry> geom, const int nstart, const int nfence, std::shared_ptr<const ZCoeff_Block> coeff,
                  const bool gaunt, const bool breit, const bool store_c = false, const bool store_g = false)
     : RelJop(geom, nstart, nfence, coeff, gaunt, breit, store_c, store_g) {
    }
};

}

#endif

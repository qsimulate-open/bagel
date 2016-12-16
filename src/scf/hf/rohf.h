//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: rohf.h
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


#ifndef __BAGEL_SRC_SCF_ROHF_H
#define __BAGEL_SRC_SCF_ROHF_H

#include <src/scf/hf/uhf.h>

// implements UHF as in Tsuchimochi and Scuseria, J. Chem. Phys. 133, 141102 (2010)
namespace bagel {

class ROHF : public UHF {
  protected:
    void symmetrize_cv(std::shared_ptr<Matrix>, std::shared_ptr<Matrix>);

  public:
    ROHF() { }
    ROHF(const std::shared_ptr<const PTree> idata, const std::shared_ptr<const Geometry> geom,
         const std::shared_ptr<const Reference> re = nullptr) : UHF(idata, geom, re) {
    }

    void compute() override;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::ROHF)

#endif

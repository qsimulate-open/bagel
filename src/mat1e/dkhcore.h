//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcore.h
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Raymond Wang <raymondwang@u.northwestern.edu>
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


#ifndef __SRC_MAT1E_DKHCORE_H
#define __SRC_MAT1E_DKHCORE_H

#include <src/util/constants.h>
#include <src/util/math/matrix.h>
#include <src/molecule/molecule.h>

namespace bagel {

class DKHcore : public Matrix {
  protected:
    std::shared_ptr<const Molecule> mol_;
    void init(std::shared_ptr<const Molecule>);

  public:
    DKHcore() { }
    DKHcore(std::shared_ptr<const Molecule>);
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::DKHcore)

#endif


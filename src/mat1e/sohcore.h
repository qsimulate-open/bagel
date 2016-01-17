//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: sohcore.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __SRC_MAT1E_SOHCORE_H
#define __SRC_MAT1E_SOHCORE_H

#include <src/util/math/zmatrix.h>
#include <src/molecule/molecule.h>
#include <src/mat1e/hcore.h>

namespace bagel {

class SOHcore : public ZMatrix {
  protected:
    std::shared_ptr<const Hcore> hcore_;

    void form_sohcore();

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<ZMatrix>(*this) & hcore_;
    }

  public:
    SOHcore() { }
    SOHcore(std::shared_ptr<const Molecule> geom, std::shared_ptr<const Hcore> h);

};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::SOHcore)

#endif


//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: sofock.h
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


#ifndef __BAGEL_SRC_SCF_SOFOCK_H
#define __BAGEL_SRC_SCF_SOFOCK_H

#include <src/df/df.h>
#include <src/wfn/geometry.h>
#include <src/integral/libint/libint.h>
#include <src/util/math/zmatrix.h>

namespace bagel {

class SOFock : public ZMatrix {
  protected:
    std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const ZMatrix> previous_;
    std::shared_ptr<const ZMatrix> socoeff_;

    void form_sofock();

    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<ZMatrix>(*this) & geom_ & previous_ & socoeff_;
    }

  public:
    SOFock(const std::shared_ptr<const Geometry> geom, const std::shared_ptr<const ZMatrix> previous, const std::shared_ptr<const ZMatrix> socoeff);

};

}
#endif

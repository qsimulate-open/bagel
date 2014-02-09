//
// BAGEL - Parallel electron correlation program.
// Filename: space.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Modified by: Michael Caldwell <caldwell@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_FCI_SPACE_H
#define __SRC_FCI_SPACE_H

#include <src/fci/space_base.h>

namespace bagel {

class HZSpace : public Space_base {
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Space_base>(*this);
    }

  public:
    HZSpace() { }
    HZSpace(const int norb, const int nelea, const int neleb, const bool compress = false, const bool mute = false);
    HZSpace(std::shared_ptr<const Determinants> det, const bool compress = false, const bool mute = false)
      : HZSpace(det->norb(), det->nelea(), det->neleb(), compress, mute) {
    }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::HZSpace)

#endif

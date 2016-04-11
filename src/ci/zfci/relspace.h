//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relspace.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Michael Caldwell <caldwell@northwestern.edu>
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


#ifndef __SRC_ZFCI_RELSPACE_H
#define __SRC_ZFCI_RELSPACE_H

#include <src/ci/fci/space_base.h>

namespace bagel {

// implements spaces that contain all determinants |PQ> for a given Kramers index -N/2 to N/2
class RelSpace : public Space_base {
  protected:
    int nele_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & boost::serialization::base_object<Space_base>(*this);
      ar & nele_;
    }

  public:
    RelSpace() { }
    RelSpace(const int norb, const int nele, const bool mute = true, const bool linkup = false);

    int nele() const { return nele_; }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::RelSpace)

#endif

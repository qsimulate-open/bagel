//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: space.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_FCI_SPACE_H
#define __SRC_FCI_SPACE_H

#include <src/ci/fci/space_base.h>

namespace bagel {

class HZSpace : public Space_base {
  protected:
    std::shared_ptr<CIStringSpace<CIStringSet<FCIString>>> spacea_;
    std::shared_ptr<CIStringSpace<CIStringSet<FCIString>>> spaceb_;

    void link();

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<Space_base>(*this) << spacea_ << spaceb_;
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<Space_base>(*this) >> spacea_ >> spaceb_;
      link();
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

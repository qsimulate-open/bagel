//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: space_base.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_FCI_SPACE_BASE_H
#define __SRC_FCI_SPACE_BASE_H

#include <src/ci/fci/determinants.h>

namespace bagel {

class Space_base {
  protected:
    std::map<std::pair<int, int>, std::shared_ptr<Determinants>> detmap_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & detmap_;
    }

  public:
    Space_base() { }
    virtual ~Space_base() { }

    // static constants
    static const int Alpha = 0;
    static const int Beta = 1;

    std::shared_ptr<Determinants> finddet(const int i, const int j) { return detmap_.at({i,j}); }
    std::shared_ptr<const Determinants> finddet(const int i, const int j) const { return detmap_.at({i,j}); }

    const std::map<std::pair<int,int>, std::shared_ptr<Determinants>>& detmap() const { return detmap_; }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Space_base)

namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<Space_base, T>::value>::type> {
    typedef Space_base type;
  };
}

#endif

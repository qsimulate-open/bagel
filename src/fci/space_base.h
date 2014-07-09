//
// BAGEL - Parallel electron correlation program.
// Filename: space_base.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Modified by: Shane Parker <shane.parker@u.northwestern.edu>
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


#ifndef __SRC_FCI_SPACE_BASE_H
#define __SRC_FCI_SPACE_BASE_H

#include <src/fci/determinants.h>

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

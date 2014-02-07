//
// BAGEL - Parallel electron correlation program.
// Filename: determinants_base.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_CIUTIL_DETERMINANTS_BASE_H
#define __SRC_CIUTIL_DETERMINANTS_BASE_H

#include <tuple>
#include <string>
#include <map>
#include <algorithm>
#include <cassert>
#include <src/util/constants.h>
#include <src/util/serialization.h>
#include <src/ciutil/cistring.h>
#include <src/ciutil/ciblock.h>

namespace bagel {

// Determinants class without linkages.
class Determinants_base : public CIBlockInfo<FCIString> {
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
    }

  public:
    Determinants_base() { }
    Determinants_base(const int norb, const int nelea, const int neleb, const bool mute = false);
    virtual ~Determinants_base() { }

    size_t ncsfs() const;

    std::pair<std::vector<std::tuple<int, int, int>>, double> spin_adapt(const int, std::bitset<nbit__>, std::bitset<nbit__>) const;
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Determinants_base)
namespace bagel {
  template <class T>
  struct base_of<T, typename std::enable_if<std::is_base_of<Determinants_base, T>::value>::type> {
    typedef Determinants_base type;
  };
}

#endif

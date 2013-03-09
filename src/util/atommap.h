//
// BAGEL - Parallel electron correlation program.
// Filename: atommap.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __SRC_UTIL_ATOMMAP_H
#define __SRC_UTIL_ATOMMAP_H

#include <map>
#include <string>

namespace bagel {

struct AtomMap {
  public:
    AtomMap();

    std::map<std::string, int> atommap;
    std::map<std::string, double> bsradii; 
    std::map<std::string, int> angmap;

    int angular_number(const std::string) const;
    int max_angular_number() const { return angmap.size()-1; };

    int atom_number(const std::string) const;
    double radius(const std::string) const;

    const std::string angular_string(const int);
};

}

#endif

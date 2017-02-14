//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: atommap.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_UTIL_ATOMMAP_H
#define __SRC_UTIL_ATOMMAP_H

#include <map>
#include <string>
#include <tuple>

namespace bagel {

struct AtomMap {
  public:
    AtomMap();

    std::map<std::string, int> atommap;
    std::map<std::string, double> bsradii;
    std::map<std::string, double> cov_radii;
    std::map<std::string, double> nuclear_exponents;
    std::map<std::string, double> averaged_masses;
    std::map<std::string, int> angmap;
    std::map<std::string, double> hfccp;
    std::map<std::string, std::tuple<int,int,int,int>> nclosed;
    std::map<std::string, std::tuple<int,int,int,int>> nopen;

    int angular_number(const std::string) const;
    int max_angular_number() const { return angmap.size()-1; }

    int atom_number(const std::string) const;
    double radius(const std::string) const;
    double cov_radius(const std::string) const;
    double nuclear_exponent(const std::string) const;
    double averaged_mass(const std::string) const;

    bool hfcc_exists(const std::string) const;
    double hfcc_pfac(const std::string) const;

    std::tuple<int,int,int,int> num_closed(const std::string) const;
    std::tuple<int,int,int,int> num_open(const std::string) const;

    const std::string angular_string(const int);
};

}

#endif

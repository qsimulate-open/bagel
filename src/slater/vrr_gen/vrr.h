//
// Newint - Parallel electron correlation program.
// Filename: vrr.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <string>
#include <cmath>

class VRR {
  private:
    // target angular momentum
    int a_, c_;

    // rank of Rys quadruture
    int rank_;

    // generating functions for each case
    const std::pair<std::string, int> vrr00 (         ) const; 
    const std::pair<std::string, int> vrrn0 (const int) const; 
    const std::pair<std::string, int> vrr0m (const int) const; 
    const std::pair<std::string, int> vrr11 (         ) const; 
    const std::pair<std::string, int> vrrn1 (const int) const;
    const std::pair<std::string, int> vrr1m (const int) const;
    const std::pair<std::string, int> vrrnm (const int, const int) const;

  public:
    VRR(const int _i, const int _j): a_(_i), c_(_j), rank_(::ceil(0.5 * (_i + _j + 2))) { }; // for slater!! 
    ~VRR() { };

    const std::string dump(const std::string) const;

};


//
// BAGEL - Parallel electron correlation program.
// Filename: shell_base.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_MOLECULE_SHELL_BASE_H
#define __SRC_MOLECULE_SHELL_BASE_H

#include <array>

namespace bagel {

class Shell_base {

  protected:
    bool spherical_;
    std::array<double,3> position_;
    int angular_number_;

  public:
    Shell_base() { }
    Shell_base(const bool spherical, const std::array<double,3>& position, int angular_num);
    Shell_base(const bool sph);
    virtual ~Shell_base() { }

    bool spherical() const { return spherical_; };

    double position(const int i) const { return position_[i]; };
    const std::array<double,3>& position() const { return position_; };

    int angular_number() const { return angular_number_; };

    virtual std::string show() const;

};

}

#endif


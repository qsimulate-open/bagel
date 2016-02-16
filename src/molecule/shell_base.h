//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: shell_base.h
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


#ifndef __SRC_MOLECULE_SHELL_BASE_H
#define __SRC_MOLECULE_SHELL_BASE_H

#include <array>
#include <src/util/serialization.h>

namespace bagel {

class Shell_base {

  protected:
    bool spherical_;
    std::array<double,3> position_;
    int angular_number_;

  private:
    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & spherical_ & position_ & angular_number_;
    }

  public:
    Shell_base() { }
  protected:
    Shell_base(const bool spherical, const std::array<double,3>& position, int angular_num);
    Shell_base(const bool sph);

  public:
    virtual ~Shell_base() { }

    bool spherical() const { return spherical_; }

    double position(const int i) const { return position_[i]; }
    const std::array<double,3>& position() const { return position_; }

    int angular_number() const { return angular_number_; }

    virtual std::string show() const;

};

}

#endif


//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: constraint.h
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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

#include <string>
#include <array>
#include <src/util/timer.h>
#include <src/util/io/moldenout.h>
#include <src/wfn/construct_method.h>

#ifndef __SRC_OPT_CONSTRAINT_H
#define __SRC_OPT_CONSTRAINT_H

namespace bagel  {
  class OptConstraint {
    protected:
      std::string type_;
      std::array<int,4> pair_;
      double value_;

    public:
      OptConstraint(std::shared_ptr<const PTree> inp);

      const std::string type() const { return type_; }
      const std::array<int,4>& pair() const { return pair_; }
      int pair(const unsigned int i) const { return pair_[i]; }
      double value() const { return value_; }
  };
}

#endif

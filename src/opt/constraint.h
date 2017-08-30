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


#ifndef __SRC_OPT_CONSTRAINT_H
#define __SRC_OPT_CONSTRAINT_H

#include <array>
#include <src/util/input/input.h>

namespace bagel  {
  class OptExpBonds {
    protected:
      std::array<int,2> pair_;

    public:
      OptExpBonds(std::shared_ptr<const PTree> inp);

      const std::array<int,2>& pair() const { return pair_; }
      int pair(const unsigned int i) const { return pair_[i]; }
  };
}

#endif

//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: integral.h
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


// Base class for all integral classes

#ifndef __SRC_INTEGRAL_INTEGRAL_H
#define __SRC_INTEGRAL_INTEGRAL_H

namespace bagel {

enum class Int_t { Standard, London };

template <typename DataType>
class Integral_base {
  protected:

  public:
    Integral_base () {}

    virtual void compute() = 0;

    virtual DataType* data(const int i) = 0;

};

using Integral = Integral_base<double>;

}

#endif


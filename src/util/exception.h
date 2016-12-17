//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: exception.h
// Copyright (C) 2016 Toru Shiozaki
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

#ifndef __SRC_UTIL_EXCEPTION_H
#define __SRC_UTIL_EXCEPTION_H

#include <string>
#include <stdexcept>

namespace bagel {

class Termination : public std::exception {
  protected:
    const std::string message_;
  public:
    Termination(const std::string& m) : message_(m) { }
    const char* what() const noexcept override { return message_.c_str(); }
};

}

#endif

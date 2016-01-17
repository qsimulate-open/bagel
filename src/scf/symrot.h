//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: symrot.h
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


#ifndef __src_scf_symrot_h
#define __src_scf_symrot_h

#include <vector>

namespace bagel {

class SymRotAbel {
  protected:
    std::vector<std::vector<double>> primrot_;

  public:
    std::vector<double> primrot(const int i) const { return primrot_[i]; }

    SymRotAbel(const std::vector<double>&, const int, const bool);
    ~SymRotAbel();


};

}

#endif


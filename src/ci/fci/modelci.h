//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: modelci.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __BAGEL_FCI_MODELCI_H
#define __BAGEL_FCI_MODELCI_H

#include <src/util/math/matrix.h>
#include <src/ci/fci/mofile.h>

namespace bagel {

class CIHamiltonian : public Matrix {
  private:
    using SD = std::pair<std::bitset<nbit__>, std::bitset<nbit__>>;

  protected:
    std::vector<SD> basis_;
    std::shared_ptr<const MOFile> jop_;

  public:
    CIHamiltonian(std::vector<SD> b, std::shared_ptr<const MOFile> jop);
};

class CISpin : public Matrix {
  private:
    using SD = std::pair<std::bitset<nbit__>, std::bitset<nbit__>>;

  protected:
    std::vector<SD> basis_;

  public:
    CISpin(std::vector<SD>& b, const int norb);
};

}

#endif

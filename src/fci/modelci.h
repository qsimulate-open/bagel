
// BAGEL - Parallel electron correlation program.
// Filename: modelci.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __BAGEL_FCI_MODELCI_H
#define __BAGEL_FCI_MODELCI_H

#include <src/math/matrix.h>
#include <src/fci/mofile.h>

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

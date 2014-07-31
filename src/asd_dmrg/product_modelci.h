
// BAGEL - Parallel electron correlation program.
// Filename: product_modelci.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __BAGEL_ASD_DMRG_PRODUCT_MODELCI_H
#define __BAGEL_ASD_DMRG_PRODUCT_MODELCI_H

#include <src/math/matrix.h>
#include <src/dimer/dimer_jop.h>
#include <src/asd_dmrg/dmrg_block.h>

namespace bagel {

namespace PCI {
  struct Basis : public BlockKey {
    std::bitset<nbit__> alpha;
    std::bitset<nbit__> beta;
    int state;

    Basis(const BlockKey K, const int i, const std::bitset<nbit__> a, const std::bitset<nbit__> b) : BlockKey(K), alpha(a), beta(b), state(i) {}
    BlockKey key() const { return BlockKey(this->nelea, this->neleb); }
  };
}


class ProductCIHamiltonian : public Matrix {
  protected:
    std::vector<PCI::Basis> basis_;
    std::shared_ptr<const DMRG_Block> left_;
    std::shared_ptr<const DimerJop> jop_;

  public:
    ProductCIHamiltonian(std::vector<PCI::Basis>& b, std::shared_ptr<const DMRG_Block> left, std::shared_ptr<const DimerJop> jop);
};

class ProductCISpin : public Matrix {
  protected:
    std::vector<PCI::Basis> basis_;
    std::shared_ptr<const DMRG_Block> left_;

  public:
    ProductCISpin(std::vector<PCI::Basis>& b, std::shared_ptr<const DMRG_Block> left, const int norb);
};

}

#endif

//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg/block_operators.h
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

#ifndef __BAGEL_ASD_DMRG_BLOCK_OPERATORS_H
#define __BAGEL_ASD_DMRG_BLOCK_OPERATORS_H

#include <src/dimer/dimer_jop.h>

#include <src/asd_dmrg/product_civec.h>
#include <src/asd_dmrg/dmrg_block.h>

namespace bagel {

/// Storage of operators already contracted with integrals, same as in Kurashige and Yanai paper (JCP 2009 130, 234114)
class BlockOperators {
  protected:
    // Quantities keyed by only a single BlockKey do not change the block (i.e., they don't change the quantum numbers)
    /// /f[ \hat H^{L'L} = \sum_{p,q} h_{pq} \hat p^\dagger \hat q + \frac{1}{2} \sum_{pqrs} (pq|rs)p^\dagger r^\dagger s q \f]
    std::map<BlockKey, std::shared_ptr<const Matrix>> ham_;

    /// \f[ \hat Q_{r_\alpha s_\alpha}^{L'L} = \sum_{p,q} \left\{ \Gamma^{L'L}_{p^\dagger_\alpha q_\alpha} \left[(pq|rs) - (ps|qr)\right]
    ///                                   + \Gamma^{L'L}_{p^\dagger_\beta q_\beta}(pq|rs)\right\} \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> Q_aa_;

    /// \f[ \hat Q_{r_\beta s_\beta}^{L'L} = \sum_{p,q} \left\{ \Gamma^{L'L}_{p^\dagger_\beta q_\beta} \left[(pq|rs) - (ps|qr)\right]
    ///                                   + \Gamma^{L'L}_{p^\dagger_\alpha q_\alpha}(pq|rs)\right\} \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> Q_bb_;

  public:
    BlockOperators(std::shared_ptr<const DMRG_Block> left, std::shared_ptr<DimerJop> jop);

    std::shared_ptr<const Matrix> ham(const BlockKey bk) const { return ham_.at(bk); }
    std::shared_ptr<const btas::Tensor4<double>> Q_aa(const BlockKey bk) const { return Q_aa_.at(bk); }
    std::shared_ptr<const btas::Tensor4<double>> Q_bb(const BlockKey bk) const { return Q_bb_.at(bk); }
};


}

#endif


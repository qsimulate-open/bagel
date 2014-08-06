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

/// Storage of operators already contracted with integrals, same as in Kurashige and Yanai paper (JCP 2009 130, 234114).
/// All quantities are keyed by the ket block
class BlockOperators {
  protected:
    /// \f[ \hat H^{L'L} = \sum_{p,q} h_{pq} \hat p^\dagger \hat q + \frac{1}{2} \sum_{pqrs} (pq|rs)p^\dagger r^\dagger s q \f]
    std::map<BlockKey, std::shared_ptr<const Matrix>> ham_;

    /// \f[ \hat Q_{r_\alpha s_\alpha}^{L'L} = \sum_{p,q} \left\{ \Gamma^{L'L}_{p^\dagger_\alpha q_\alpha} \left[(pq|rs) - (ps|qr)\right] + \Gamma^{L'L}_{p^\dagger_\beta q_\beta}(pq|rs)\right\} \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> Q_aa_;

    /// \f[ \hat Q_{r_\beta s_\beta}^{L'L} = \sum_{p,q} \left\{ \Gamma^{L'L}_{p^\dagger_\beta q_\beta} \left[(pq|rs) - (ps|qr)\right] + \Gamma^{L'L}_{p^\dagger_\alpha q_\alpha}(pq|rs)\right\} \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> Q_bb_;

    /// \f[ \hat Q_{r_\alpha s_\beta}^{L'L} = - \sum_{p,q} \Gamma^{L'L}_{p^\dagger_\beta q_\alpha} (ps|qr) \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> Q_ab_;

    /// \f[ \hat S_{p_\alpha}^{L'L} = \sum_q h_{pq} \Gamma^{L'L}_{q^\dagger_\alpha} + \sum_{\tau}\sum_{q,r,s} \Gamma^{L'L}_{s^\dagger_\alpha r^\dagger_\tau q_\tau} (ps|rq)  \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor3<double>>> S_a_;

    /// \f[ \hat S_{p_\beta}^{L'L} = \sum_q h_{pq} \Gamma^{L'L}_{q^\dagger_\beta} + \sum_{\tau}\sum_{q,r,s} \Gamma^{L'L}_{s^\dagger_\beta r^\dagger_\tau q_\tau} (ps|rq) \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor3<double>>> S_b_;

    /// \f[ \hat P_{r_\alpha s_\alpha}^{L'L} = \sum_{p,q} \Gamma^{L'L}_{p_\alpha q_\alpha} (pr|qs) \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> P_aa_;

    /// \f[ \hat P_{r_\beta s_\beta}^{L'L} = \sum_{p,q} \Gamma^{L'L}_{p_\beta q_\beta} (pr|qs) \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> P_bb_;

    /// \f[ \hat P_{r_\alpha s_\beta}^{L'L} = \sum_{p,q} \Gamma^{L'L}_{p_\beta q_\alpha} (pr|qs) \f]
    std::map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> P_ab_;

    /// \f[ \hat D_{p_\alpha q_\tau q_\tau}^{L'L} = \sum_s \Gamma^{L'L}_{s_\alpha} (ps|rq) \f]

    /** This is memory inefficient but convenient for now */
    std::map<BlockKey, std::shared_ptr<const btas::TensorN<double,5>>> D_a_;

    /// \f[ \hat D_{p_\beta q_\tau r_\tau}^{L'L} = \sum_s \Gamma^{L'L}_{s_\beta} (ps|rq) \f]

    /** This is memory inefficient but convenient for now */
    std::map<BlockKey, std::shared_ptr<const btas::TensorN<double,5>>> D_b_;

  public:
    BlockOperators(std::shared_ptr<const DMRG_Block> left, std::shared_ptr<DimerJop> jop);

    std::shared_ptr<const Matrix> ham(const BlockKey bk) const { return ham_.at(bk); }
    std::shared_ptr<const btas::Tensor4<double>> Q_aa(const BlockKey bk) const { return Q_aa_.at(bk); }
    std::shared_ptr<const btas::Tensor4<double>> Q_bb(const BlockKey bk) const { return Q_bb_.at(bk); }
    std::shared_ptr<const btas::Tensor4<double>> Q_ab(const BlockKey bk) const { return Q_ab_.at(bk); }
    std::shared_ptr<const btas::Tensor3<double>> S_a(const BlockKey bk) const { return S_a_.at(bk); }
    std::shared_ptr<const btas::Tensor3<double>> S_b(const BlockKey bk) const { return S_b_.at(bk); }
    std::shared_ptr<const btas::TensorN<double,5>> D_a(const BlockKey bk) const { return D_a_.at(bk); }
    std::shared_ptr<const btas::TensorN<double,5>> D_b(const BlockKey bk) const { return D_b_.at(bk); }
};

}

#endif


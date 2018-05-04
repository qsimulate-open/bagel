//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg/block_operators.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __BAGEL_ASD_DMRG_BLOCK_OPERATORS_H
#define __BAGEL_ASD_DMRG_BLOCK_OPERATORS_H

#include <src/asd/dimer/dimer_jop.h>
#include <src/asd/dmrg/block_key.h>
#include <src/util/math/blocksparsematrix.h>

namespace bagel {

class DMRG_Block;
class DMRG_Block1;
class DMRG_Block2;

/// Base class for access to contracted operators. All operators are keyed by the ket of the operation. The contracted operators that get built are the following:
/** \f[ \hat H^{L'L} = \sum_{p,q} h_{pq} \hat p^\dagger \hat q + \frac{1}{2} \sum_{pqrs} (pq|rs)p^\dagger r^\dagger s q \f]
    \f[ \hat Q_{r_\alpha^\dagger s_\alpha}^{L'L} = \sum_{p,q} \left\{ \Gamma^{L'L}_{p^\dagger_\alpha q_\alpha} \left[(pq|rs) - (ps|qr)\right] + \Gamma^{L'L}_{p^\dagger_\beta q_\beta}(pq|rs)\right\} \f]
    \f[ \hat Q_{r_\beta^\dagger s_\beta}^{L'L} = \sum_{p,q} \left\{ \Gamma^{L'L}_{p^\dagger_\beta q_\beta} \left[(pq|rs) - (ps|qr)\right] + \Gamma^{L'L}_{p^\dagger_\alpha q_\alpha}(pq|rs)\right\} \f]
    \f[ \hat Q_{r_\alpha^\dagger s_\beta}^{L'L} = - \sum_{p,q} \Gamma^{L'L}_{p^\dagger_\beta q_\alpha} (ps|qr) \f]
    \f[ \hat S_{p_\alpha}^{L'L} = \sum_q h_{pq} \Gamma^{L'L}_{q^\dagger_\alpha} + \sum_{\tau}\sum_{q,r,s} \Gamma^{L'L}_{s^\dagger_\alpha r^\dagger_\tau q_\tau} (ps|rq)  \f]
    \f[ \hat S_{p_\beta}^{L'L} = \sum_q h_{pq} \Gamma^{L'L}_{q^\dagger_\beta} + \sum_{\tau}\sum_{q,r,s} \Gamma^{L'L}_{s^\dagger_\beta r^\dagger_\tau q_\tau} (ps|rq) \f]
    \f[ \hat P_{r_\alpha^\dagger s_\alpha^\dagger}^{L'L} = \sum_{p,q} \Gamma^{L'L}_{p_\alpha q_\alpha} (pr|qs) \f]
    \f[ \hat P_{r_\beta^\dagger s_\beta^\dagger}^{L'L} = \sum_{p,q} \Gamma^{L'L}_{p_\beta q_\beta} (pr|qs) \f]
    \f[ \hat P_{r_\alpha^\dagger s_\beta^\dagger}^{L'L} = \sum_{p,q} \Gamma^{L'L}_{p_\beta q_\alpha} (pr|qs) \f]
    \f[ \hat D_{q_\tau^\dagger r_\tau s_\alpha}^{L'L} = \sum_s \Gamma^{L'L}_{p_\alpha^\dagger} (ps|rq) \f]
    \f[ \hat D_{q_\tau^\dagger r_\tau s_\beta}^{L'L} = \sum_s \Gamma^{L'L}_{p_\beta^\dagger} (ps|rq) \f] */
class BlockOperators {
  protected:
    double thresh_; ///< threshold for building BlockSparseMatrix objects

  public:
    BlockOperators(const double thresh) : thresh_(thresh) {}

    virtual std::shared_ptr<Matrix>     ham(const BlockKey bk) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix>    Q_aa(const BlockKey bk, int i, int j) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix>    Q_bb(const BlockKey bk, int i, int j) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix>    Q_ab(const BlockKey bk, int i, int j) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix>     S_a(const BlockKey bk, int i) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix>     S_b(const BlockKey bk, int i) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix>    P_aa(const BlockKey bk, int i, int j) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix>    P_bb(const BlockKey bk, int i, int j) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix>    P_ab(const BlockKey bk, int i, int j) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix> gamma_a(const BlockKey bk, int i) const  = 0;
    virtual std::shared_ptr<BlockSparseMatrix> gamma_b(const BlockKey bk, int i) const  = 0;

    virtual double  ham(const BlockKey bk, int brastate, int ketstate) const  = 0;
    virtual double Q_aa(const BlockKey bk, int brastate, int ketstate, int i, int j) const  = 0;
    virtual double Q_bb(const BlockKey bk, int brastate, int ketstate, int i, int j) const  = 0;
    virtual double Q_ab(const BlockKey bk, int brastate, int ketstate, int i, int j) const  = 0;
    virtual double  S_a(const BlockKey bk, int brastate, int ketstate, int i) const  = 0;
    virtual double  S_b(const BlockKey bk, int brastate, int ketstate, int i) const  = 0;
    virtual double P_aa(const BlockKey bk, int brastate, int ketstate, int i, int j) const  = 0;
    virtual double P_bb(const BlockKey bk, int brastate, int ketstate, int i, int j) const  = 0;
    virtual double P_ab(const BlockKey bk, int brastate, int ketstate, int i, int j) const  = 0;
    virtual double  D_a(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const  = 0;
    virtual double  D_b(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const  = 0;
};

/// Operators for a single block
class BlockOperators1 : public BlockOperators {
  protected:
    std::shared_ptr<const DMRG_Block1> left_;

    std::unordered_map<BlockKey, std::shared_ptr<const Matrix>> ham_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> Q_aa_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> Q_bb_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> Q_ab_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::Tensor3<double>>> S_a_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::Tensor3<double>>> S_b_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> P_aa_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> P_bb_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::Tensor4<double>>> P_ab_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::TensorN<double,5>>> D_a_;
    std::unordered_map<BlockKey, std::shared_ptr<const btas::TensorN<double,5>>> D_b_;

  public:
    BlockOperators1(std::shared_ptr<const DMRG_Block1> left, std::shared_ptr<DimerJop> jop, const double thresh = 1.0e-12);

    std::shared_ptr<Matrix>     ham(const BlockKey bk) const override { return ham_.at(bk)->copy(); }
    std::shared_ptr<BlockSparseMatrix>    Q_aa(const BlockKey bk, int i, int j) const override { return get_sparse_mat_block(Q_aa_.at(bk), i, j); }
    std::shared_ptr<BlockSparseMatrix>    Q_bb(const BlockKey bk, int i, int j) const override { return get_sparse_mat_block(Q_bb_.at(bk), i, j); }
    std::shared_ptr<BlockSparseMatrix>    Q_ab(const BlockKey bk, int i, int j) const override { return get_sparse_mat_block(Q_ab_.at(bk), i, j); }
    std::shared_ptr<BlockSparseMatrix>     S_a(const BlockKey bk, int i) const override { return get_sparse_mat_block(S_a_.at(bk), i); }
    std::shared_ptr<BlockSparseMatrix>     S_b(const BlockKey bk, int i) const override { return get_sparse_mat_block(S_b_.at(bk), i); }
    std::shared_ptr<BlockSparseMatrix>    P_aa(const BlockKey bk, int i, int j) const override { return get_sparse_mat_block(P_aa_.at(bk), i, j); }
    std::shared_ptr<BlockSparseMatrix>    P_bb(const BlockKey bk, int i, int j) const override { return get_sparse_mat_block(P_bb_.at(bk), i, j); }
    std::shared_ptr<BlockSparseMatrix>    P_ab(const BlockKey bk, int i, int j) const override { return get_sparse_mat_block(P_ab_.at(bk), i, j); }
    std::shared_ptr<BlockSparseMatrix> gamma_a(const BlockKey bk, int i) const override;
    std::shared_ptr<BlockSparseMatrix> gamma_b(const BlockKey bk, int i) const override;

    std::shared_ptr<Matrix> gamma_a_as_matrix(const BlockKey bk, int i) const;
    std::shared_ptr<Matrix> gamma_b_as_matrix(const BlockKey bk, int i) const;

    void  S_a_copy_to(double* target, const BlockKey bk, int i) const { copy_to(target, S_a_.at(bk), i); }
    void  S_b_copy_to(double* target, const BlockKey bk, int i) const { copy_to(target, S_b_.at(bk), i); }
    void Q_aa_copy_to(double* target, const BlockKey bk, int i, int j) const { copy_to(target, Q_aa_.at(bk), i, j); }
    void Q_bb_copy_to(double* target, const BlockKey bk, int i, int j) const { copy_to(target, Q_bb_.at(bk), i, j); }
    void Q_ab_copy_to(double* target, const BlockKey bk, int i, int j) const { copy_to(target, Q_ab_.at(bk), i, j); }
    void P_aa_copy_to(double* target, const BlockKey bk, int i, int j) const { copy_to(target, P_aa_.at(bk), i, j); }
    void P_bb_copy_to(double* target, const BlockKey bk, int i, int j) const { copy_to(target, P_bb_.at(bk), i, j); }
    void P_ab_copy_to(double* target, const BlockKey bk, int i, int j) const { copy_to(target, P_ab_.at(bk), i, j); }

    double  ham(const BlockKey bk, int brastate, int ketstate) const override { return ham_.at(bk)->element(brastate,ketstate); }
    double Q_aa(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return (*Q_aa_.at(bk))(brastate, ketstate, i, j); }
    double Q_bb(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return (*Q_bb_.at(bk))(brastate, ketstate, i, j); }
    double Q_ab(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return (*Q_ab_.at(bk))(brastate, ketstate, i, j); }
    double  S_a(const BlockKey bk, int brastate, int ketstate, int i) const override { return (*S_a_.at(bk))(brastate, ketstate, i); }
    double  S_b(const BlockKey bk, int brastate, int ketstate, int i) const override { return (*S_b_.at(bk))(brastate, ketstate, i); }
    double P_aa(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return (*P_aa_.at(bk))(brastate, ketstate, i, j); }
    double P_bb(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return (*P_bb_.at(bk))(brastate, ketstate, i, j); }
    double P_ab(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return (*P_ab_.at(bk))(brastate, ketstate, i, j); }
    double  D_a(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const override { return (*D_a_.at(bk))(brastate, ketstate, i, j, k); }
    double  D_b(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const override { return (*D_b_.at(bk))(brastate, ketstate, i, j, k); }

    const MatView Q_aa_as_matview(const BlockKey bk) const { return get_as_matview(Q_aa_.at(bk)); }
    const MatView Q_bb_as_matview(const BlockKey bk) const { return get_as_matview(Q_bb_.at(bk)); }
    const MatView Q_ab_as_matview(const BlockKey bk) const { return get_as_matview(Q_ab_.at(bk)); }
    const MatView  S_a_as_matview(const BlockKey bk) const { return get_as_matview(S_a_.at(bk)); }
    const MatView  S_b_as_matview(const BlockKey bk) const { return get_as_matview(S_b_.at(bk)); }
    const MatView P_aa_as_matview(const BlockKey bk) const { return get_as_matview(P_aa_.at(bk)); }
    const MatView P_bb_as_matview(const BlockKey bk) const { return get_as_matview(P_bb_.at(bk)); }
    const MatView P_ab_as_matview(const BlockKey bk) const { return get_as_matview(P_ab_.at(bk)); }
    const MatView  D_a_as_matview(const BlockKey bk) const { return get_as_matview(D_a_.at(bk)); }
    const MatView  D_b_as_matview(const BlockKey bk) const { return get_as_matview(D_b_.at(bk)); }

  private:
    template <typename TensorType>
    const MatView get_as_matview(std::shared_ptr<const TensorType> t) const {
      return MatView(btas::make_view(btas::CRange<2>(t->extent(0)*t->extent(1), t->size()/(t->extent(0)*t->extent(1))), t->storage()), true);
    }

    template <typename... Args>
    std::shared_ptr<BlockSparseMatrix> get_sparse_mat_block(Args... args) const {
      return std::make_shared<BlockSparseMatrix>(get_mat_block(args...));
    }

    std::shared_ptr<Matrix> get_mat_block(std::shared_ptr<const btas::Tensor3<double>> g, const int i) const {
      auto out = std::make_shared<Matrix>(g->extent(0), g->extent(1), true);
      std::copy_n(&(*g)(0, 0, i), out->size(), out->data());
      return out;
    }

    std::shared_ptr<Matrix> get_mat_block(std::shared_ptr<const btas::Tensor4<double>> g, const int i, const int j) const {
      auto out = std::make_shared<Matrix>(g->extent(0), g->extent(1), true);
      std::copy_n(&(*g)(0, 0, i, j), out->size(), out->data());
      return out;
    }

    std::shared_ptr<Matrix> get_mat_block(std::shared_ptr<const btas::TensorN<double,5>> g, const int i, const int j, const int k) const {
      auto out = std::make_shared<Matrix>(g->extent(0), g->extent(1), true);
      std::copy_n(&(*g)(0, 0, i, j, k), out->size(), out->data());
      return out;
    }

    void copy_to(double* target, std::shared_ptr<const btas::Tensor3<double>> g, const int i) const {
      std::copy_n(&(*g)(0, 0, i), g->extent(0)*g->extent(1), target);
    }

    void copy_to(double* target, std::shared_ptr<const btas::Tensor4<double>> g, const int i, const int j) const {
      std::copy_n(&(*g)(0, 0, i, j), g->extent(0)*g->extent(1), target);
    }

};

/// Operators for two blocks
class BlockOperators2 : public BlockOperators {
  protected:
    std::shared_ptr<const DMRG_Block2> blocks_;
    std::shared_ptr<const BlockOperators1> left_ops_;
    std::shared_ptr<const BlockOperators1> right_ops_;
    std::shared_ptr<const BlockOperators1> intra_ops_;
    std::shared_ptr<DimerJop> jop_;

  public:
    BlockOperators2(std::shared_ptr<const DMRG_Block2> blocks, std::shared_ptr<DimerJop> jop, const double thresh = 1.0e-13);

    std::shared_ptr<Matrix>     ham(const BlockKey bk) const override;
    std::shared_ptr<BlockSparseMatrix>    Q_aa(const BlockKey bk, int i, int j) const override;
    std::shared_ptr<BlockSparseMatrix>    Q_bb(const BlockKey bk, int i, int j) const override;
    std::shared_ptr<BlockSparseMatrix>    Q_ab(const BlockKey bk, int i, int j) const override;
    std::shared_ptr<BlockSparseMatrix>     S_a(const BlockKey bk, int i) const override;
    std::shared_ptr<BlockSparseMatrix>     S_b(const BlockKey bk, int i) const override;
    std::shared_ptr<BlockSparseMatrix>    P_aa(const BlockKey bk, int i, int j) const override;
    std::shared_ptr<BlockSparseMatrix>    P_bb(const BlockKey bk, int i, int j) const override;
    std::shared_ptr<BlockSparseMatrix>    P_ab(const BlockKey bk, int i, int j) const override;
    std::shared_ptr<BlockSparseMatrix> gamma_a(const BlockKey bk, int i) const override;
    std::shared_ptr<BlockSparseMatrix> gamma_b(const BlockKey bk, int i) const override;

    // these should all be replaced with more direct functions
    double  ham(const BlockKey bk, int brastate, int ketstate) const override { return ham(bk)->element(brastate,ketstate); }
    double Q_aa(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return Q_aa(bk, i, j)->element(brastate, ketstate); }
    double Q_bb(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return Q_bb(bk, i, j)->element(brastate, ketstate); }
    double Q_ab(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return Q_ab(bk, i, j)->element(brastate, ketstate); }
    double  S_a(const BlockKey bk, int brastate, int ketstate, int i) const override { return S_a(bk, i)->element(brastate, ketstate); }
    double  S_b(const BlockKey bk, int brastate, int ketstate, int i) const override { return S_b(bk, i)->element(brastate, ketstate); }
    double P_aa(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return P_aa(bk, i, j)->element(brastate, ketstate); }
    double P_bb(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return P_bb(bk, i, j)->element(brastate, ketstate); }
    double P_ab(const BlockKey bk, int brastate, int ketstate, int i, int j) const override { return P_ab(bk, i, j)->element(brastate, ketstate); }
    double  D_a(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const override;
    double  D_b(const BlockKey bk, int brastate, int ketstate, int i, int j, int k) const override;
};

}

#endif


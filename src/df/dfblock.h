//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dfblock.h
// Copyright (C) 2012 Toru Shiozaki
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

#ifndef __SRC_DF_DFBLOCK_H
#define __SRC_DF_DFBLOCK_H

#include <numeric>
#include <src/util/math/algo.h>
#include <src/util/timer.h>
#include <src/util/simple.h>
#include <src/util/taskqueue.h>
#include <src/util/parallel/staticdist.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/math/btas_interface.h>
#include <src/util/math/matrix.h>
#include <src/util/math/matop.h>

namespace bagel {

/*
    DFBlock is a slice of 3-index DF integrals. Distributed by the first index
*/

class DFBlock : public btas::Tensor3<double> {

  // aux_ runs fastest, b2_ runs slowest
  public:
    using btas::Tensor3<double>::data;

  protected:
    // distribution information
    std::shared_ptr<const StaticDist> adist_shell_;
    std::shared_ptr<const StaticDist> adist_;

    // if true, asize is evenly distributed. If false, asize is at the shell boundary
    bool averaged_;

    // a set of offsets of this block in the entire DF integrals
    size_t astart_;
    size_t b1start_;
    size_t b2start_;

  public:

    DFBlock() { }

    // construction of a block from AO integrals
    DFBlock(std::shared_ptr<const StaticDist> adist_shell, std::shared_ptr<const StaticDist> adist,
                 const size_t a, const size_t b1, const size_t b2, const int as, const int b1s, const int b2s, const bool averaged = false);
    DFBlock(const DFBlock& o);

    // dimensions of the block
    size_t asize() const { return this->extent(0); }
    size_t b1size() const { return this->extent(1); }
    size_t b2size() const { return this->extent(2); }

    size_t size() const { return asize()*b1size()*b2size(); }
    bool averaged() const { return averaged_; }

    // a set of offsets of this block in the entire DF integrals
    size_t astart() const { return astart_; }
    size_t b1start() const { return b1start_; }
    size_t b2start() const { return b2start_; }

    // dist
    const std::shared_ptr<const StaticDist>& adist_now() const { return averaged_ ? adist_ : adist_shell_; }


    // some math functions
    DFBlock& operator=(const DFBlock& o);
    DFBlock& operator=(DFBlock&& o);
    DFBlock& operator+=(const DFBlock& o) { btas::Tensor3<double>::operator+=(o); return *this; }
    DFBlock& operator-=(const DFBlock& o) { btas::Tensor3<double>::operator-=(o); return *this; }

    template <typename ScaleType, class DType>
    void ax_plus_y(const ScaleType a, const DType& o) { btas::axpy(a, o, *this); }
    template <typename ScaleType, class DType>
    void ax_plus_y(const ScaleType a, const std::shared_ptr<DType>& o) { ax_plus_y(a, *o); }
    template <typename ScaleType>
    void scale(const ScaleType a) { btas::scal(a, *this); }

    void zero() { std::fill_n(data(), size(), 0.0); }

    // symmetrize b1 and b2 (assuming b1size() == b2size())
    void symmetrize();

    void copy_block(std::shared_ptr<MatView> o, const int jdim, const size_t offset);
    void copy_block(MatView o, const int jdim, const size_t offset);
    void add_block(std::shared_ptr<MatView> o, const int jdim, const size_t offset, const double fac = 1.0);
    void add_block(MatView o, const int jdim, const size_t offset, const double fac = 1.0);

    // average the asize between MPI processes (block will be described by dist_)
    void average();

    // reverse operation of average() function
    void shell_boundary();

    std::shared_ptr<DFBlock> clone() const;
    std::shared_ptr<DFBlock> copy() const;

    std::shared_ptr<DFBlock> transform_second(const MatView c, const bool trans = false) const;
    std::shared_ptr<DFBlock> transform_third(const MatView c, const bool trans = false) const;

    std::shared_ptr<DFBlock> merge_b1(std::shared_ptr<const DFBlock> o) const;
    std::shared_ptr<DFBlock> slice_b1(const int start, const int size) const;

    // add ab^+  to this.
    void add_direct_product(const VecView a, const MatView b, const double fac);

    // exchange b1 and b2
    std::shared_ptr<DFBlock> swap() const;

    // 2RDM contractions
    std::shared_ptr<DFBlock> apply_rhf_2RDM(const double scale_exch) const;
    std::shared_ptr<DFBlock> apply_uhf_2RDM(const btas::Tensor2<double>&, const btas::Tensor2<double>&) const;
    std::shared_ptr<DFBlock> apply_2RDM(const btas::Tensor4<double>& rdm, const btas::Tensor2<double>& rdm1, const int nclosed, const int nact) const;
    std::shared_ptr<DFBlock> apply_2RDM(const btas::Tensor4<double>& rdm) const;
    std::shared_ptr<DFBlock> apply_2RDM_tr(const btas::Tensor4<double>& rdm, const btas::Tensor2<double>& rdm1, const int nclosed, const int nact) const;

    // Form 2- and 4-index integrals
    std::shared_ptr<Matrix> form_2index(const std::shared_ptr<const DFBlock> o, const double a) const;
    std::shared_ptr<Matrix> form_4index(const std::shared_ptr<const DFBlock> o, const double a) const;
    // slowest index of o is fixed to n
    std::shared_ptr<Matrix> form_4index_1fixed(const std::shared_ptr<const DFBlock> o, const double a, const size_t n) const;
    std::shared_ptr<Matrix> form_4index_diagonal() const;
    std::shared_ptr<Matrix> form_4index_diagonal_part() const;
    std::shared_ptr<Matrix> form_aux_2index(const std::shared_ptr<const DFBlock> o, const double a) const;

    std::shared_ptr<VectorB> form_vec(const std::shared_ptr<const Matrix> den) const;
    std::shared_ptr<Matrix> form_mat(const btas::Tensor1<double>& fit) const;

    void contrib_apply_J(const std::shared_ptr<const DFBlock> o, const std::shared_ptr<const Matrix> mat);

    // compute (D|ia)(ia|j) and set to the location specified by the offset
    std::shared_ptr<Matrix> form_Dj(const std::shared_ptr<const Matrix> o, const int jdim) const;

    // CAUTION, ist, jst, and kst are absolute number (NOT relative to astart_, ...). Returns double[] whose size is i*j*k
    std::shared_ptr<btas::Tensor3<double>> get_block(const int ist, const int i, const int jst, const int j, const int kst, const int k) const;

};

}


#endif

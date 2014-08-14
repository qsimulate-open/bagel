//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef __SRC_DF_DFBLOCK_H
#define __SRC_DF_DFBLOCK_H

#include <numeric>
#include <src/math/algo.h>
#include <src/util/timer.h>
#include <src/util/simple.h>
#include <src/util/taskqueue.h>
#include <src/parallel/staticdist.h>
#include <src/parallel/mpi_interface.h>
#include <src/math/btas_interface.h>
#include <src/math/matrix.h>
#include <src/math/matop.h>

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
    DFBlock& operator=(const DFBlock& o) {
      btas::Tensor3<double>::operator=(o);
      adist_shell_ = o.adist_shell_;
      adist_ = o.adist_;
      averaged_ = o.averaged_;
      astart_ = o.astart_;
      b1start_ = o.b1start_;
      b2start_ = o.b2start_;
      return *this;
    }
    DFBlock& operator=(DFBlock&& o) {
      btas::Tensor3<double>::operator=(std::move(o));
      adist_shell_ = o.adist_shell_;
      adist_ = o.adist_;
      averaged_ = o.averaged_;
      astart_ = o.astart_;
      b1start_ = o.b1start_;
      b2start_ = o.b2start_;
      return *this;
    }
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
    void symmetrize() {
      if (b1size() != b2size()) throw std::logic_error("illegal call of DFBlock::symmetrize()");
      const int n = b1size();
      for (int i = 0; i != n; ++i)
        for (int j = i; j != n; ++j) {
          blas::ax_plus_y_n(1.0, data()+asize()*(j+n*i), asize(), data()+asize()*(i+n*j));
          std::copy_n(data()+asize()*(i+n*j), asize(), data()+asize()*(j+n*i));
        }
    }

    void copy_block(std::shared_ptr<MatView> o, const int jdim, const size_t offset) {
      assert(o->size() == asize()*jdim);
      std::copy_n(o->data(), asize()*jdim, data()+offset);
    }

    void copy_block(MatView o, const int jdim, const size_t offset) {
      assert(o.size() == asize()*jdim);
      std::copy_n(o.data(), asize()*jdim, data()+offset);
    }

    void add_block(std::shared_ptr<MatView> o, const int jdim, const size_t offset, const double fac = 1.0) {
      assert(o->size() == asize()*jdim);
      blas::ax_plus_y_n(fac, o->data(), asize()*jdim, data()+offset);
    }

    void add_block(MatView o, const int jdim, const size_t offset, const double fac = 1.0) {
      assert(o.size() == asize()*jdim);
      blas::ax_plus_y_n(fac, o.data(), asize()*jdim, data()+offset);
    }

    // average the asize between MPI processes (block will be described by dist_)
    void average();

    // reverse operation of average() function
    void shell_boundary();

    std::shared_ptr<DFBlock> clone() const;
    std::shared_ptr<DFBlock> copy() const;

    std::shared_ptr<DFBlock> transform_second(const MatView c, const bool trans = false) const;
    std::shared_ptr<DFBlock> transform_third(const MatView c, const bool trans = false) const;

    // add ab^+  to this.
    void add_direct_product(const std::shared_ptr<const VectorB> a, const std::shared_ptr<const Matrix> b, const double fac);

    // exchange b1 and b2
    std::shared_ptr<DFBlock> swap() const;

    // 2RDM contractions
    std::shared_ptr<DFBlock> apply_rhf_2RDM(const double scale_exch) const;
    std::shared_ptr<DFBlock> apply_uhf_2RDM(const btas::Tensor2<double>&, const btas::Tensor2<double>&) const;
    std::shared_ptr<DFBlock> apply_2RDM(const btas::Tensor4<double>& rdm, const btas::Tensor2<double>& rdm1, const int nclosed, const int nact) const;
    std::shared_ptr<DFBlock> apply_2RDM(const btas::Tensor4<double>& rdm) const;

    // Form 2- and 4-index integrals
    std::shared_ptr<Matrix> form_2index(const std::shared_ptr<const DFBlock> o, const double a) const;
    std::shared_ptr<Matrix> form_4index(const std::shared_ptr<const DFBlock> o, const double a) const;
    // slowest index of o is fixed to n
    std::shared_ptr<Matrix> form_4index_1fixed(const std::shared_ptr<const DFBlock> o, const double a, const size_t n) const;
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

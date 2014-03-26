//
// BAGEL - Parallel electron correlation program.
// Filename: paralleldf.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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


#ifndef __SRC_DF_PARALLELDF_H
#define __SRC_DF_PARALLELDF_H

#include <src/df/dfinttask_old.h>
#include <src/df/dfinttask.h>
#include <src/df/dfblock.h>
#include <src/df/dfblock_london.h>
#include <src/math/matrix.h>
#include <src/math/zmatrix.h>

namespace bagel {

template <typename DataType, typename MatrixType, typename BlockType>
class ParallelDFit : public std::enable_shared_from_this<ParallelDFit<DataType,MatrixType,BlockType>> {
  protected:
    // blocks that this process has
    std::vector<std::shared_ptr<BlockType>> block_;

    // naux runs fastest, nindex2 runs slowest
    const size_t naux_;
    const size_t nindex1_;
    const size_t nindex2_;

    std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> df_;
    // data2_ is usually empty (except for the original DFDist)
    // AO two-index integrals ^ -1/2
    std::shared_ptr<MatrixType> data2_;

    bool serial_;

  public:
    ParallelDFit<DataType,MatrixType,BlockType>(const size_t, const size_t, const size_t, std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> = nullptr, std::shared_ptr<MatrixType> = nullptr);

    size_t naux() const { return naux_; }
    size_t nindex1() const { return nindex1_; }
    size_t nindex2() const { return nindex2_; }
    size_t size() const { return naux_*nindex1_*nindex2_; }

    std::vector<std::shared_ptr<BlockType>>& block() { return block_; }
    const std::vector<std::shared_ptr<BlockType>>& block() const { return block_; }
    std::shared_ptr<BlockType> block(const size_t i) { return block_[i]; }
    std::shared_ptr<const BlockType> block(const size_t i) const { return block_[i]; }

    std::shared_ptr<const StaticDist> adist_now() const { return block_[0]->adist_now(); }

    void add_block(std::shared_ptr<BlockType> o);

    std::shared_ptr<MatrixType> form_2index(std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o, const double a, const bool swap = false) const;
    std::shared_ptr<MatrixType> form_4index(std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o, const double a, const bool swap = false) const;
    std::shared_ptr<MatrixType> form_aux_2index(std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o, const double a) const;

    void ax_plus_y(const double a, const std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o);
    void scale(const double a);

    std::shared_ptr<MatrixType> get_block(const int i, const int id, const int j, const int jd, const int k, const int kd) const;

    const std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> df() const { return df_; }
    std::shared_ptr<const MatrixType> data2() const { return data2_; }

    // compute a J operator, given density matrices in AO basis
    std::shared_ptr<MatrixType> compute_Jop(const std::shared_ptr<const MatrixType> den) const;
    std::shared_ptr<MatrixType> compute_Jop(const std::shared_ptr<const ParallelDFit<DataType,MatrixType,BlockType>> o, const std::shared_ptr<const MatrixType> den, const bool onlyonce = false) const;
    std::shared_ptr<MatrixType> compute_Jop_from_cd(std::shared_ptr<const MatrixType> cd) const;
    std::shared_ptr<MatrixType> compute_cd(const std::shared_ptr<const MatrixType> den, std::shared_ptr<const MatrixType> dat2 = nullptr, const bool onlyonce = false) const;

    void average_3index() {
      Timer time;
      if (!serial_)
        for (auto& i : block_)
          i->average();
      time.tick_print("3-index ints post");
    }

    void shell_boundary_3index() {
      if (!serial_)
        for (auto& i : block_)
          i->shell_boundary();
    }
};

using ParallelDF = ParallelDFit<double, Matrix, DFBlock>;

}

#define PARALLELDF_HEADERS
#include <src/df/paralleldf_impl.hpp>
#undef PARALLELDF_HEADERS

extern template class bagel::ParallelDFit<double, bagel::Matrix, bagel::DFBlock>;
extern template class bagel::ParallelDFit<std::complex<double>, bagel::ZMatrix, bagel::DFBlock_London>;

#endif


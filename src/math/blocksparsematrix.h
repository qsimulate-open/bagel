//
// BAGEL - Parallel electron correlation program.
// Filename: blocksparsematrix.h
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __math_blocksparsematrix_h
#define __math_blocksparsematrix_h

#include <src/math/matrix.h>

/// Class for blockwise sparse matrices
namespace bagel {

class BlockSparseMatrix {
  protected:
    /// Store blocks as separate Matrix objects with offsets
    std::map<std::pair<size_t, size_t>, std::shared_ptr<Matrix>> data_;

    int ndim_;
    int mdim_;

  public:
    /// A single block (so actually a dense matrix)
    BlockSparseMatrix(std::shared_ptr<Matrix> m);
    /// Pre-determined blocking structure
    BlockSparseMatrix(const int n, const int m, std::map<std::pair<size_t, size_t>, std::shared_ptr<Matrix>> d);
    // for completeness, a constructor that accepts a matrix and a threshold should probably exist at some point

    std::map<std::pair<size_t, size_t>, std::shared_ptr<Matrix>>& data() { return data_; }
    const std::map<std::pair<size_t, size_t>, std::shared_ptr<Matrix>>& data() const { return data_; }

    /// get block based on starting location. this is the best way to alter the data
    std::shared_ptr<Matrix> get_block(size_t nstart, size_t mstart) { return (data_.find({nstart, mstart})!=data_.end() ? data_[{nstart,mstart}] : nullptr); }
    /// const get block based on starting location.
    std::shared_ptr<const Matrix> get_block(size_t nstart, size_t mstart) const { return (data_.find({nstart, mstart})!=data_.end() ? data_.at({nstart,mstart}) : nullptr); }

    int ndim() const { return ndim_; }
    int mdim() const { return mdim_; }
    int size() const { return accumulate(data_.begin(), data_.end(), 0, [] (int x, std::pair<std::pair<size_t, size_t>, std::shared_ptr<Matrix>> p) { return x + p.second->size(); }); }

    /// element-wise read-only: is slow and not recommended.
    double element(const int n, const int m) const {
      auto iter = std::find_if(data_.begin(), data_.end(), [&n, &m] (std::pair<std::pair<size_t, size_t>, std::shared_ptr<Matrix>> p)
        { return (n >= p.first.first && n < p.first.first + p.second->ndim()) && (m >= p.first.second && m < p.first.second + p.second->mdim()); });
      if (iter!=data_.end())
        return iter->second->element(n - iter->first.first, m - iter->first.second);
      else
        return 0.0;
    }

    /// returns a VectorB containing the diagonal elements
    std::shared_ptr<VectorB> diagonal() const;
};

void mat_block_multiply(bool Atrans, bool Btrans, double alpha, const Matrix& A, const BlockSparseMatrix& B, double beta, Matrix& C);

}

#endif

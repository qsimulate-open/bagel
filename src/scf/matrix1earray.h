//
// BAGEL - Parallel electron correlation program.
// Filename: matrix1earray.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __src_scf_matrix1earray_h
#define __src_scf_matrix1earray_h

#include <cassert>
#include <src/wfn/shell.h>
#include <src/wfn/geometry.h>
#include <src/scf/matrix1e.h>
#include <string>
#include <algorithm>
#include <memory>

namespace bagel {

// specialized matrix for N component 1e integrals
template <int N>
class Matrix1eArray : public Matrix1e {
  protected:
    std::array<std::shared_ptr<Matrix>, N> matrices_;

    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int) = 0;

  public:
    Matrix1eArray(const std::shared_ptr<const Geometry>);
    Matrix1eArray(const std::shared_ptr<const Geometry>, const int n, const int m);
    Matrix1eArray(const Matrix1eArray&);

    const std::shared_ptr<const Geometry> geom() const { return geom_; }

    std::shared_ptr<Matrix>& data(const int i) { return matrices_[i]; }
    constexpr const int nblocks() { return N; }

    void fill_upper() override { for (int i = 0 ; i < N; ++i) matrices_[i]->fill_upper(); }

};

template <int N>
Matrix1eArray<N>::Matrix1eArray(const std::shared_ptr<const Geometry> geom) : Matrix1e(geom) {
  for(int i = 0; i < N; ++i) {
    matrices_[i] = std::shared_ptr<Matrix>(new Matrix(geom->nbasis(), geom->nbasis()));
    matrices_[i]->zero();
  }
}


template <int N>
Matrix1eArray<N>::Matrix1eArray(const Matrix1eArray& o) : Matrix1eArray(o.geom()) {
  for (int i = 0; i < N; ++i) {
    copy_n(o.data(i)->data(), ndim_*mdim_, this->data(i)->data()); 
  }
}

}

#endif

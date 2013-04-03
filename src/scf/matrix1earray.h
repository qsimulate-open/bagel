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
class Matrix1eArray {
  protected:
    std::array<std::shared_ptr<Matrix>, N> matrices_;

    std::shared_ptr<const Geometry> geom_;

    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int) = 0;

  public:
    Matrix1eArray(const std::shared_ptr<const Geometry>);
    Matrix1eArray(const std::shared_ptr<const Geometry>, const int n, const int m);
    Matrix1eArray(const Matrix1eArray&);

    const std::shared_ptr<const Geometry> geom() const { return geom_; }

    std::shared_ptr<Matrix>& data(const int i) { return matrices_[i]; }
    std::shared_ptr<Matrix>& operator[](const int i) { return data(i); }
    std::shared_ptr<const Matrix> data(const int i) const { return matrices_[i]; }
    std::shared_ptr<const Matrix> operator[](const int i) const { return matrices_[i]; }
    constexpr const int nblocks() { return N; }

    void fill_upper() { for (int i = 0 ; i < N; ++i) matrices_[i]->fill_upper(); }

    virtual void print(const std::string name = "") const;
    virtual void init();

};

template <int N>
Matrix1eArray<N>::Matrix1eArray(const std::shared_ptr<const Geometry> geom) : geom_(geom) {
  for(int i = 0; i < N; ++i) {
    matrices_[i] = std::shared_ptr<Matrix>(new Matrix(geom->nbasis(), geom->nbasis()));
  }
}


template <int N>
Matrix1eArray<N>::Matrix1eArray(const std::shared_ptr<const Geometry> geom, const int n, const int m) : geom_(geom) {
  for(int i = 0; i < N; ++i) {
    matrices_[i] = std::shared_ptr<Matrix>(new Matrix(n, m));
  }
}


template <int N>
Matrix1eArray<N>::Matrix1eArray(const Matrix1eArray& o) : Matrix1eArray(o.geom()) {
  const int ndim = matrices_.front()->ndim();
  const int mdim = matrices_.front()->mdim();
  for (int i = 0; i < N; ++i) {
    copy_n(o.data(i)->data(), ndim*mdim, this->data(i)->data()); 
  }
}


template <int N>
void Matrix1eArray<N>::print(const std::string name) const {
  int j = 0;
  for (auto& i : matrices_) {
    std::stringstream ss; ss << name << " " << j++;
    i->print(ss.str());
  }
}

template <int N>
void Matrix1eArray<N>::init() {

  // identical to Matrix1e::init()
  // only lower half will be stored
  // TODO rewrite. thread. parallel. distribute

  auto o0 = geom_->offsets().begin();
  int u = 0;
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++o0) {
    // iatom1 = iatom1;
    auto offset0 = o0->begin();
    for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++offset0) {
      auto offset1 = o0->begin();
      for (auto b1 = (*a0)->shells().begin(); b1 != (*a0)->shells().end(); ++b1, ++offset1) {
        if (u++ % mpi__->size() == mpi__->rank()) {
          std::array<std::shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          computebatch(input, *offset0, *offset1);
        }   
      }   
    }   

    auto o1 = o0+1;
    for (auto a1 = a0+1; a1 != geom_->atoms().end(); ++a1, ++o1) {
      auto offset0 = o0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++offset0) {
        auto offset1 = o1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++offset1) {
          if (u++ % mpi__->size() == mpi__->rank()) {
            std::array<std::shared_ptr<const Shell>,2> input = {{*b1, *b0}};
            computebatch(input, *offset0, *offset1);
          }   
        }   
      }   
    }   
  }
  //mpi__->allreduce(data_.get(), size());

}

}

#endif

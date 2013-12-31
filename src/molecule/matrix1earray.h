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


#ifndef __SRC_MOLECULE_MATRIX1EARRAY_H
#define __SRC_MOLECULE_MATRIX1EARRAY_H

#include <src/molecule/matrix1e.h>

namespace bagel {

// specialized matrix for N component 1e integrals
template <int N>
class Matrix1eArray {
  protected:
    std::array<std::shared_ptr<Matrix>, N> matrices_;

    std::shared_ptr<const Molecule> mol_;

    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int) = 0;

    bool localized_;

  public:
    Matrix1eArray(const std::shared_ptr<const Molecule>, const bool loc = false);
    Matrix1eArray(const std::shared_ptr<const Molecule>, const int n, const int m, const bool loc = false);
    Matrix1eArray(const Matrix1eArray&);

    void ax_plus_y(const double a, const Matrix1eArray<N>& o) {
      std::transform(o.matrices_.begin(), o.matrices_.end(), matrices_.begin(), matrices_.begin(),
                     [&a](std::shared_ptr<Matrix> p, std::shared_ptr<Matrix> q) { q->ax_plus_y(a, p); return q; });
    }

    const std::shared_ptr<const Molecule> mol() const { return mol_; }

    std::shared_ptr<Matrix>& data(const int i) { return matrices_[i]; }
    std::shared_ptr<Matrix>& operator[](const int i) { return data(i); }
    std::shared_ptr<const Matrix> data(const int i) const { return matrices_[i]; }
    std::shared_ptr<const Matrix> operator[](const int i) const { return matrices_[i]; }
    constexpr static int Nblocks() { return N; }

    void fill_upper() { for (int i = 0 ; i < N; ++i) matrices_[i]->fill_upper(); }

    virtual void print(const std::string name = "") const;
    virtual void init();

    void localize() {
      localized_ = true;
      for (auto& i : matrices_) i->localize();
    }

};

template <int N>
Matrix1eArray<N>::Matrix1eArray(const std::shared_ptr<const Molecule> mol, const bool loc) : mol_(mol), localized_(loc) {
  static_assert(N > 0, "Matrix1eArray should be constructed with N > 0");
  for(int i = 0; i < N; ++i) {
    matrices_[i] = std::make_shared<Matrix>(mol->nbasis(), mol->nbasis());
  }
}


template <int N>
Matrix1eArray<N>::Matrix1eArray(const std::shared_ptr<const Molecule> mol, const int n, const int m, const bool loc) : mol_(mol), localized_(loc) {
  static_assert(N > 0, "Matrix1eArray should be constructed with N > 0");
  for(int i = 0; i < N; ++i) {
    matrices_[i] = std::make_shared<Matrix>(n, m);
  }
}


template <int N>
Matrix1eArray<N>::Matrix1eArray(const Matrix1eArray& o) : Matrix1eArray(o.mol(), o.localized_) {
  for (int i = 0; i < N; ++i) {
    *data(i) = *o.data(i);
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

  size_t oa0 = 0;
  int u = 0;
  for (auto a0 = mol_->atoms().begin(); a0 != mol_->atoms().end(); ++a0) {
    // iatom1 = iatom1;
    size_t ob0 = oa0;
    for (auto& b0 : (*a0)->shells()) {
      size_t ob1 = oa0;
      for (auto& b1 : (*a0)->shells()) {
        if (u++ % mpi__->size() == mpi__->rank()) {
          computebatch({{b1, b0}}, ob0, ob1);
        }
        ob1 += b1->nbasis();
      }
      ob0 += b0->nbasis();
    }

    auto oa1 = oa0 + (*a0)->nbasis();
    for (auto a1 = a0+1; a1 != mol_->atoms().end(); ++a1) {
      size_t ob0 = oa0;
      for (auto& b0 : (*a0)->shells()) {
        size_t ob1 = oa1;
        for (auto& b1 : (*a1)->shells()) {
          if (u++ % mpi__->size() == mpi__->rank()) {
            computebatch({{b1, b0}}, ob0, ob1);
          }
          ob1 += b1->nbasis();
        }
        ob0 += b0->nbasis();
      }
      oa1 += (*a1)->nbasis();
    }
    oa0 += (*a0)->nbasis();
  }
  for (auto& i : matrices_) i->allreduce();

}

}

#endif

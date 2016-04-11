//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: matrix1earray.h
// Copyright (C) 2013 Toru Shiozaki
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


#ifndef __SRC_MOLECULE_MATRIX1EARRAY_H
#define __SRC_MOLECULE_MATRIX1EARRAY_H

#include <src/mat1e/matrix1e.h>
#include <src/util/taskqueue.h>

namespace bagel {

template <int N, typename MatType>
class Matrix1eArrayTask;

// specialized matrix for N component 1e integrals
template <int N, typename MatType = Matrix>
class Matrix1eArray {
  friend class Matrix1eArrayTask<N, MatType>;
  protected:
    std::array<std::shared_ptr<MatType>, N> matrices_;

    virtual void init(std::shared_ptr<const Molecule>);
    virtual void computebatch(const std::array<std::shared_ptr<const Shell>,2>&, const int, const int, std::shared_ptr<const Molecule>) = 0;

    bool localized_;

  public:
    Matrix1eArray() {}
    Matrix1eArray(const std::shared_ptr<const Molecule>, const bool loc = false);
    Matrix1eArray(const int n, const int m, const bool loc = false);
    Matrix1eArray(const Matrix1eArray&);
    virtual ~Matrix1eArray() { }

    void ax_plus_y(const double a, const Matrix1eArray<N, MatType>& o) {
      std::transform(o.matrices_.begin(), o.matrices_.end(), matrices_.begin(), matrices_.begin(),
                     [&a](std::shared_ptr<MatType> p, std::shared_ptr<MatType> q) { q->ax_plus_y(a, p); return q; });
    }

    std::shared_ptr<MatType>& data(const int i) { return matrices_[i]; }
    std::shared_ptr<const MatType> data(const int i) const { return matrices_[i]; }

    MatType& operator[](const int i) { return *matrices_[i]; }
    const MatType& operator[](const int i) const { return *matrices_[i]; }

    constexpr static int Nblocks() { return N; }

    void fill_upper()       { for (auto& i : matrices_) i->fill_upper(); }
    void fill_upper_conjg() { for (auto& i : matrices_) i->fill_upper_conjg(); }
    template<typename DataType>
    void scale(const DataType a) { for (auto& i : matrices_) i->scale(a); }

    virtual void print(const std::string name = "", const int len = 10) const;

    void localize() {
      localized_ = true;
      for (auto& i : matrices_) i->localize();
    }

};

template <int N, typename MatType>
Matrix1eArray<N, MatType>::Matrix1eArray(const std::shared_ptr<const Molecule> mol, const bool loc) : localized_(loc) {
  static_assert(N > 0, "Matrix1eArray should be constructed with N > 0");
  for(int i = 0; i < N; ++i) {
    matrices_[i] = std::make_shared<MatType>(mol->nbasis(), mol->nbasis(), loc);
  }
}


template <int N, typename MatType>
Matrix1eArray<N, MatType>::Matrix1eArray(const int n, const int m, const bool loc) : localized_(loc) {
  static_assert(N > 0, "Matrix1eArray should be constructed with N > 0");
  for(int i = 0; i < N; ++i) {
    matrices_[i] = std::make_shared<MatType>(n, m, loc);
  }
}


template <int N, typename MatType>
Matrix1eArray<N, MatType>::Matrix1eArray(const Matrix1eArray& o) : localized_(o.localized_) {
  for (int i = 0; i < N; ++i) {
    *data(i) = *o.data(i);
  }
}


template <int N, typename MatType>
void Matrix1eArray<N, MatType>::print(const std::string name, const int len) const {
  int j = 0;
  for (auto& i : matrices_) {
    std::stringstream ss; ss << name << " " << j++;
    i->print(ss.str(), len);
  }
}

template <int N, typename MatType>
void Matrix1eArray<N, MatType>::init(std::shared_ptr<const Molecule> mol) {

  // identical to Matrix1e::init()
  // only lower half will be stored
  const size_t nshell = accumulate(mol->atoms().begin(), mol->atoms().end(), 0, [](int r, std::shared_ptr<const Atom> p) { return r+p->nshell(); });
  TaskQueue<Matrix1eArrayTask<N, MatType>> task(nshell*(nshell+1)/2);

  size_t oa0 = 0;
  int u = 0;
  for (auto a0 = mol->atoms().begin(); a0 != mol->atoms().end(); ++a0) {
    // iatom1 = iatom1;
    size_t ob0 = oa0;
    for (auto& b0 : (*a0)->shells()) {
      size_t ob1 = oa0;
      for (auto& b1 : (*a0)->shells()) {
        if (u++ % mpi__->size() == mpi__->rank()) {
          task.emplace_back(std::array<std::shared_ptr<const Shell>,2>{{b1, b0}}, ob0, ob1, mol, this);
        }
        ob1 += b1->nbasis();
      }
      ob0 += b0->nbasis();
    }

    auto oa1 = oa0 + (*a0)->nbasis();
    for (auto a1 = a0+1; a1 != mol->atoms().end(); ++a1) {
      size_t ob0 = oa0;
      for (auto& b0 : (*a0)->shells()) {
        size_t ob1 = oa1;
        for (auto& b1 : (*a1)->shells()) {
          if (u++ % mpi__->size() == mpi__->rank()) {
            task.emplace_back(std::array<std::shared_ptr<const Shell>,2>{{b1, b0}}, ob0, ob1, mol, this);
          }
          ob1 += b1->nbasis();
        }
        ob0 += b0->nbasis();
      }
      oa1 += (*a1)->nbasis();
    }
    oa0 += (*a0)->nbasis();
  }
  task.compute();
  for (auto& i : matrices_) i->allreduce();

}

template <int N, typename MatType = Matrix>
class Matrix1eArrayTask {
  protected:
    Matrix1eArray<N, MatType>* parent_;
    size_t ob0, ob1;
    std::array<std::shared_ptr<const Shell>,2> bas;
    std::shared_ptr<const Molecule> mol;
  public:
    Matrix1eArrayTask<N, MatType>(std::array<std::shared_ptr<const Shell>,2> a, size_t b, size_t c, std::shared_ptr<const Molecule> m, Matrix1eArray<N, MatType>* d)
      : parent_(d), ob0(b), ob1(c), bas(a), mol(m) { }
    void compute() const { parent_->computebatch(bas, ob0, ob1, mol); }
};

}

#endif

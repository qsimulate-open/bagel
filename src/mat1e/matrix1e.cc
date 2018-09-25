//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: matrix1e.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <src/mat1e/matrix1e.h>
#include <src/util/f77.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/parallel/resources.h>
#include <src/util/taskqueue.h>

using namespace std;
using namespace bagel;

template <typename MatType, class Enable>
Matrix1e_<MatType, Enable>::Matrix1e_(shared_ptr<const Molecule> mol) : MatType(mol->nbasis(), mol->nbasis()) {
  zero();
}


template <typename MatType, class Enable>
Matrix1e_<MatType, Enable>::Matrix1e_(const Matrix1e_<MatType, Enable>& o) : MatType(o.ndim(), o.mdim()) {
  copy_n(o.data(), size(), data());
}


namespace bagel {
template <typename MatType, class Enable = typename std::enable_if<(std::is_same<MatType, Matrix>::value || std::is_same<MatType, ZMatrix>::value)>::type>
class Matrix1eTask_ {
  protected:
    Matrix1e_<MatType>* parent_;
    size_t ob0, ob1;
    array<shared_ptr<const Shell>,2> bas;
    shared_ptr<const Molecule> mol;
  public:
    Matrix1eTask_(array<shared_ptr<const Shell>,2> a, size_t b, size_t c, shared_ptr<const Molecule> m, Matrix1e_<MatType>* d)
      : parent_(d), ob0(b), ob1(c), bas(a), mol(m) { }
    void compute() const { parent_->computebatch(bas, ob0, ob1, mol); }
};
}


template <typename MatType, class Enable>
void Matrix1e_<MatType, Enable>::init(shared_ptr<const Molecule> mol) {

  // CAUTION only lower half will be stored
  const size_t nshell = accumulate(mol->atoms().begin(), mol->atoms().end(), 0, [](int r, shared_ptr<const Atom> p) { return r+p->nshell(); });
  TaskQueue<Matrix1eTask_<MatType, Enable>> task(nshell*(nshell+1)/2);

  size_t oa0 = 0;
  int u = 0;
  for (auto a0 = mol->atoms().begin(); a0 != mol->atoms().end(); ++a0) {
    // iatom1 = iatom1;
    size_t ob0 = oa0;
    for (auto& b0 : (*a0)->shells()) {
      size_t ob1 = oa0;
      for (auto& b1 : (*a0)->shells()) {
        if (u++ % mpi__->size() == mpi__->rank()) {
          task.emplace_back(array<shared_ptr<const Shell>,2>{{b1, b0}}, ob0, ob1, mol, this);
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
            task.emplace_back(array<shared_ptr<const Shell>,2>{{b1, b0}}, ob0, ob1, mol, this);
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
  allreduce();
}


template class bagel::Matrix1e_<Matrix>;
template class bagel::Matrix1e_<ZMatrix>;


BOOST_CLASS_EXPORT_IMPLEMENT(Matrix1e)
BOOST_CLASS_EXPORT_IMPLEMENT(ZMatrix1e)



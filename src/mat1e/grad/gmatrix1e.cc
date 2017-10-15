//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gmatrix1e.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#include <src/mat1e/grad/gmatrix1e.h>
#include <src/util/taskqueue.h>

using namespace std;
using namespace bagel;

template <typename MatType>
GMatrix1e<MatType>::GMatrix1e(const shared_ptr<const Molecule> mol) {
  matrices_.assign(3 * mol->natom(), make_shared<MatType>(mol->nbasis(), mol->nbasis()));
}


template <typename MatType>
GMatrix1e<MatType>::GMatrix1e(const GMatrix1e& o) {
  matrices_.resize(o.size());
  for (int i = 0; i < size(); ++i) {
    *data(i) = *o.data(i);
  }
}


template <typename MatType>
void GMatrix1e<MatType>::print(const string name, const int len) const {
  int j = 0;
  for (auto& i : matrices_) {
    stringstream ss; ss << name << " " << j++;
    i->print(ss.str(), len);
  }
}


template <typename MatType>
void GMatrix1e<MatType>::init(shared_ptr<const Molecule> mol) {

  vector<shared_ptr<GMatrix1eTask<Matrix>>> task;
  const size_t nshell  = accumulate(mol->atoms().begin(), mol->atoms().end(), 0,
                                          [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
  task.reserve(nshell*nshell);

  // TODO perhaps we could reduce operation by a factor of 2
  int cnt = 0;
  int iatom0 = 0;
  auto oa0 = mol->offsets().begin();
  for (auto a0 = mol->atoms().begin(); a0 != mol->atoms().end(); ++a0, ++oa0, ++iatom0) {
    int iatom1 = 0;
    auto oa1 = mol->offsets().begin();
    for (auto a1 = mol->atoms().begin(); a1 != mol->atoms().end(); ++a1, ++oa1, ++iatom1) {

      auto o0 = oa0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++o0) {
        auto o1 = oa1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++o1) {

          // static distribution since this is cheap
          if (cnt++ % mpi__->size() != mpi__->rank()) continue;

          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          vector<int> atom = {iatom0, iatom1};
          vector<int> offset = {*o0, *o1};

          task.push_back(make_shared<GMatrix1eTask<Matrix>>(input, atom, offset, mol, this));
        }
      }
    }
  }

  TaskQueue<shared_ptr<GMatrix1eTask<Matrix>>> tq(move(task));
  tq.compute();

  for (auto& i : matrices_) i->allreduce();

}


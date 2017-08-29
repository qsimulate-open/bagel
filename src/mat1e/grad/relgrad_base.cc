//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relgrad_base.cc
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


#include <src/mat1e/grad/relgrad_base.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/parallel/resources.h>
#include <src/util/math/xyzfile.h>
#include <src/util/taskqueue.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(Relgrad_base)


Relgrad_base::Relgrad_base(shared_ptr<const Molecule> mol0)
  : Matrix(mol0->nbasis(), mol0->nbasis()), nbasis(mol0->nbasis()) mol(mol0)) {

  molu = make_shared<Molecule>(*mol);
  molu = molu->uncontract();
  nunc = molu->nbasis();
  U = make_shared<MixedBasis<OverlapBatch>>(molu, mol);
  U_T = make_shared<MixedBasis<OverlapBatch>>(mol, molu);

  s = make_shared<Overlap>(molu);

  // id = make_shared<Matrix>(nunc, nunc);
  // for (int i = 0; i < nunc; i++) {
  //   (*id)(i, i) = 1;
  // }

  void store_mat(shared_ptr<VectorB> vec) {
    shared_ptr<Matrix> mat = make_shared<Matrix>(nunc, nunc);
    for (int i = 0; i < nunc; i++) {
      (*mat)(i, i) = *vec(i);
    }
    vec2mat[vec] = mat;
  }
}

shared_ptr<std::vector<std::vector<const Matrix>>> Relgrad_base::overlapgrad() {
  shared_ptr<vector<mutex>> mut = make_shared<vector<mutex>>(mol->natom());
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      if (k == l) {

      }
      else {
        shared_ptr<GradFile> grad = make_shared<GradFile>(mol->natom());
        grad->zero();

        shared_ptr<Matrix> den = make_shared<Matrix>(nbasis, nbasis);
        for (int m = 0; m < nbasis; k++) {
          for (int n = 0; n < nbasis; l++) {
            (*den)(m, n) = (*U_T)(k, m) * (*U)(n, l) / ((*s)(l, l) - (*s)(k, k));
          }
        }
        vector<shared_ptr<GradTask>> task;
        const size_t nshell  = std::accumulate(mol->atoms().begin(), mol->atoms().end(), 0,
                              [](const int& i, const shared_ptr<const Atom>& o) { return i+o->shells().size(); });
        task.reserve(nshell*nshell);
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
                if (cnt++ % mpi__->size() != mpi__->rank()) continue;
                array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
                vector<int> atom = {iatom0, iatom1};
                vector<int> offset_ = {*o0, *o1};
                task.push_back(make_shared<GradTask1s2>(input, atom, offset_, den, mut, grad, mol));
              }
            }
          }
        }

        TaskQueue<shared_ptr<GradTask>> tq(move(task));
        tq.compute();
        grad->allreduce();
      }
    }
  }
}

shared_ptr<std::vector<std::vector<const Matrix>>> Relgrad_base::kineticgrad() {

}

shared_ptr<std::vector<std::vector<const Matrix>>> Relgrad_base::naigrad() {

}

shared_ptr<std::vector<std::vector<const Matrix>>> Relgrad_base::smallnaigrad() {

}
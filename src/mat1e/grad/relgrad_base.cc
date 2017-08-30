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
#include <src/grad/gradeval_base.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(Relgrad_base)


Relgrad_base::Relgrad_base(shared_ptr<const Geometry> geom)
  : Matrix(geom->nbasis(), geom->nbasis()), natom(geom->natom()) nbasis(geom->nbasis()), mol(geom) {

  molu = make_shared<Molecule>(*mol->uncontract());
  nunc = molu->nbasis();
  U_T = make_shared<MixedBasis<OverlapBatch>>(mol, molu);

  s = make_shared<Overlap>(molu);

  id = make_shared<Matrix>(nunc, nunc);
  for (int i = 0; i < nunc; i++) {
    (*id)(i, i) = 1;
  }

  void store_mat(shared_ptr<VectorB> vec) {
    shared_ptr<Matrix> mat = make_shared<Matrix>(nunc, nunc);
    for (int i = 0; i < nunc; i++) {
      (*mat)(i, i) = *vec(i);
    }
    vec2mat[vec] = mat;
  }
}

shared_ptr<vector<vector<Matrix>>> Relgrad_base::overlapgrad() {
  shared_ptr<const Matrix> U = make_shared<Matrix>(U_T % id);
  shared_ptr<vector<vector<Matrix>>> PU = make_shared<vector<vector<Matrix>>>(3, vector<Matrix>(natom, Matrix(nunc, nunc)));
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      if (k == l) {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < natom; j++) {
            (*PU)[i][j](k, l) = 0;
          }
        }
      }
      else {
        shared_ptr<GradFile> grad = make_shared<GradFile>(natom);
        shared_ptr<Matrix> den = make_shared<Matrix>(nbasis, nbasis);
        for (int m = 0; m < nbasis; k++) {
          for (int n = 0; n < nbasis; l++) {
            (*den)(m, n) = (*U_T)(k, m) * (*U)(n, l) / ((*s)(l, l) - (*s)(k, k));
          }
        }
        GradEval_base ge(mol);
        grad = ge.contract_overlapgrad(den);
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < natom; j++) {
            (*PU)[i][j](k, l) = grad(i, j);
          }
        }
      }
    }
  }
}

shared_ptr<std::vector<std::vector<const Matrix>>> Relgrad_base::kineticgrad() {
  shared_ptr<Matrix> K_X = make_shared<Matrix>(Matrix(1024, 1024));
  return K_X;
}

shared_ptr<std::vector<std::vector<const Matrix>>> Relgrad_base::naigrad() {
  shared_ptr<Matrix> V_X = make_shared<Matrix>(Matrix(1024, 1024));
  return V_X;
}

shared_ptr<std::vector<std::vector<const Matrix>>> Relgrad_base::smallnaigrad() {
  shared_ptr<Matrix> O_X = make_shared<Matrix>(Matrix(1024, 1024));
  return O_X;
}
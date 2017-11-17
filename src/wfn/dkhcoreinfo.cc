//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcoreinfo.cc
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

#include <src/wfn/dkhcoreinfo.h>
#include <src/mat1e/kinetic.h>
#include <src/mat1e/nai.h>
#include <src/mat1e/rel/small1e.h>
#include <src/mat1e/overlap.h>
#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/mixedbasis.h>

using namespace std;
using namespace bagel;

DKHcoreInfo::DKHcoreInfo(const shared_ptr<const Molecule> current, const int dkh_level)
  : dkh2(false), tgrad_(make_shared<const GKinetic>(current)), vgrad_(nullptr), natom_(current->natom()), nbasis_(current->nbasis()) {
  switch (dkh_level) {
  case 2:
    vgrad_ = make_shared<const GNAI>(current);
    dkh2 = true;
    init_t(current);
    init_v(current);
    init_v2(current);
    break;
  case 1:
    vgrad_ = make_shared<const GNAI>(current);
    init_v(current);
  case 0:
    init_t(current);
    break;
  default:
    throw runtime_error("Only order 0-2 allowed for DKH.");
  }

  auto mol = make_shared<Molecule>(*current);
  mol = mol->uncontract();
  const MixedBasis<OverlapBatch> mix(mol, current);

  VectorB eig;
  const Overlap overlap(mol);
  shared_ptr<const Matrix> tildex = overlap.tildex();
  const Kinetic kinetic(mol);
  auto tmp = make_shared<Matrix>(*tildex % kinetic * *tildex);
  int nunc = tmp->ndim();
  eig = VectorB(nunc);
  tmp->diagonalize(eig);
  transfer_ = make_shared<Matrix>(mix * *tildex * *tmp);

  // check if conditions are satisfied
  const Overlap s(current);
  const Matrix s_p = *transfer_ % s * *transfer_;
  const Kinetic t(current);
  const Matrix t_p = *transfer_ % t * *transfer_;
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      if (k == l) {
        assert(s_p(k, l) == 1);
        assert(t_p(k, l) != 0);
      }
      else {
        assert(s_p(k, l) == 0);
        assert(t_p(k, l) == 0);
      }
    }
  }
  
}

void DKHcoreInfo::init_t(const shared_ptr<const Molecule> current) {
  trelgrad_ = vector<Matrix>(3 * natom_, Matrix(nbasis_, nbasis_));
}

void DKHcoreInfo::init_v(const shared_ptr<const Molecule> current) {
  vrelgrad_ = vector<Matrix>(3 * natom_, Matrix(nbasis_, nbasis_));
}

void DKHcoreInfo::init_v2(const shared_ptr<const Molecule> current) {
  v2relgrad_ = vector<Matrix>(3 * natom_, Matrix(nbasis_, nbasis_));
}

shared_ptr<GradFile> DKHcoreInfo::compute_t(const array<shared_ptr<const Shell>,2>& s, const array<int,4>& a, const array<int,4>& o, const shared_ptr<const Matrix> den) const {
  const int dimb1 = s[0]->nbasis();
  const int dimb0 = s[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(o[1], o[0], dimb1, dimb0);
  
  auto out = make_shared<GradFile>(natom);

  for (int k = 0; k < 3; k++) {
    std::shared_ptr<const Matrix> ct0 = trelgrad_[3 * a[0] + k]->get_submatrix(o[1], o[0], dimb1, dimb0);
    std::shared_ptr<const Matrix> ct1 = trelgrad_[3 * a[1] + k]->get_submatrix(o[1], o[0], dimb1, dimb0);
    out->element(k, jatom0) += ddot_(cden->size(), cden->data(), 1, ct0->data(), 1);
    out->element(k, jatom1) += ddot_(cden->size(), cden->data(), 1, ct1->data(), 1);
  }
  return out;
}

shared_ptr<GradFile> DKHcoreInfo::compute_v(const array<shared_ptr<const Shell>,2>& s, const array<int,4>& a, const array<int,4>& o, const shared_ptr<const Matrix> den) const {
  const int dimb1 = s[0]->nbasis();
  const int dimb0 = s[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(o[1], o[0], dimb1, dimb0);
  
  auto out = make_shared<GradFile>(natom);

  for (int k = 0; k < 3; k++) {
    std::shared_ptr<const Matrix> cv0 = vrelgrad_[3 * a[0] + k]->get_submatrix(o[1], o[0], dimb1, dimb0);
    std::shared_ptr<const Matrix> cv1 = vrelgrad_[3 * a[1] + k]->get_submatrix(o[1], o[0], dimb1, dimb0);
    out->element(k, jatom0) += ddot_(cden->size(), cden->data(), 1, cv0->data(), 1);
    out->element(k, jatom1) += ddot_(cden->size(), cden->data(), 1, cv1->data(), 1);
  }
  return out;
}

shared_ptr<GradFile> DKHcoreInfo::compute_v2(const array<shared_ptr<const Shell>,2>& s, const array<int,4>& a, const array<int,4>& o, const shared_ptr<const Matrix> den) const {
  const int dimb1 = s[0]->nbasis();
  const int dimb0 = s[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(o[1], o[0], dimb1, dimb0);
  
  auto out = make_shared<GradFile>(natom);

  for (int k = 0; k < 3; k++) {
    std::shared_ptr<const Matrix> cv20 = v2relgrad_[3 * a[0] + k]->get_submatrix(o[1], o[0], dimb1, dimb0);
    std::shared_ptr<const Matrix> cv21 = v2relgrad_[3 * a[1] + k]->get_submatrix(o[1], o[0], dimb1, dimb0);
    out->element(k, jatom0) += ddot_(cden->size(), cden->data(), 1, cv20->data(), 1);
    out->element(k, jatom1) += ddot_(cden->size(), cden->data(), 1, cv21->data(), 1);
  }
  return out;
}


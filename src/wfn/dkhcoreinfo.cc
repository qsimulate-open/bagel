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
  : dkh2_(false), vgrad_(nullptr), natom_(current->natom()), nbasis_(current->nbasis()) {
  mol_ = current->uncontract();
  mix_ = make_shared<MixedBasis<OverlapBatch>>(current, mol_);

  VectorB eig;
  const Overlap overlap(mol_);
  shared_ptr<const Matrix> tildex = overlap.tildex();
  const Kinetic kinetic(mol_);
  auto tmp = make_shared<Matrix>(*tildex % kinetic * *tildex);
  int nunc = tmp->ndim();
  eig = VectorB(nunc);
  tmp->diagonalize(eig);
  transfer_ = make_shared<Matrix>(*tildex * *tmp);

  // check if conditions are satisfied
  const Matrix s_p = *transfer_ % overlap * *transfer_;
  const Matrix t_p = *transfer_ % kinetic * *transfer_;
  s_p.print("s_p");
  t_p.print("t_p");
  for (int k = 0; k < nunc; k++) {
    for (int l = 0; l < nunc; l++) {
      if (k == l) {
        assert(fabs(s_p(k, l) - 1) < 1.0e-8);
        assert(fabs(t_p(k, l)) > 1.0e-8);
      }
      else {
        assert(fabs(s_p(k, l)) < 1.0e-8);
        assert(fabs(t_p(k, l)) < 1.0e-8);
      }
    }
  }

  tgrad_ = make_shared<const GKinetic>(mol_);
  switch (dkh_level) {
  case 2:
    vgrad_ = make_shared<const GNAI>(mol_);
    dkh2_ = true;
    init_t();
    init_v();
    init_v2();
    break;
  case 1:
    vgrad_ = make_shared<const GNAI>(mol_);
    init_v();
  case 0:
    init_t();
    break;
  default:
    throw runtime_error("Only order 0-2 allowed for DKH.");
  }
  
}

// Free particle FW term
void DKHcoreInfo::init_t() {
  cout << "  Analytical gradient will involve DKH0 corrections to the 1e Hamiltonian." << endl;

  trelgrad_ = vector<Matrix>(3 * natom_, Matrix(nbasis_, nbasis_));
  for (int i = 0; i < 3 * natom_; i++) {
    trelgrad_[i] = *mix_ % (*tgrad_)[i] * *mix_;
  }
}

// First order FW transformation of nuclear potential
void DKHcoreInfo::init_v() {
  cout << "  Analytical gradient will involve DKH1 corrections to the 1e Hamiltonian." << endl;

  vrelgrad_ = vector<Matrix>(3 * natom_, Matrix(nbasis_, nbasis_));
}

// Second order contributions
void DKHcoreInfo::init_v2() {
  cout << "  Analytical gradient will involve DKH2 corrections to the 1e Hamiltonian." << endl;

  v2relgrad_ = vector<Matrix>(3 * natom_, Matrix(nbasis_, nbasis_));
}

shared_ptr<GradFile> DKHcoreInfo::compute_t(const array<shared_ptr<const Shell>,2>& s, const array<int,4>& a, const array<int,4>& o, const shared_ptr<const Matrix> den) const {
  const int dimb1 = s[0]->nbasis();
  const int dimb0 = s[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(o[1], o[0], dimb1, dimb0);
  
  auto out = make_shared<GradFile>(natom_);

  for (int k = 0; k < 3; k++) {
    std::shared_ptr<const Matrix> ct0 = trelgrad_[3 * a[0] + k].get_submatrix(o[1], o[0], dimb1, dimb0);
    std::shared_ptr<const Matrix> ct1 = trelgrad_[3 * a[1] + k].get_submatrix(o[1], o[0], dimb1, dimb0);
    out->element(k, a[0]) += ddot_(cden->size(), cden->data(), 1, ct0->data(), 1);
    out->element(k, a[1]) += ddot_(cden->size(), cden->data(), 1, ct1->data(), 1);
  }
  return out;
}

shared_ptr<GradFile> DKHcoreInfo::compute_v(const array<shared_ptr<const Shell>,2>& s, const array<int,4>& a, const array<int,4>& o, const shared_ptr<const Matrix> den) const {
  assert(vgrad());
  const int dimb1 = s[0]->nbasis();
  const int dimb0 = s[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(o[1], o[0], dimb1, dimb0);
  
  auto out = make_shared<GradFile>(natom_);

  for (int k = 0; k < 3; k++) {
    std::shared_ptr<const Matrix> cv0 = vrelgrad_[3 * a[0] + k].get_submatrix(o[1], o[0], dimb1, dimb0);
    std::shared_ptr<const Matrix> cv1 = vrelgrad_[3 * a[1] + k].get_submatrix(o[1], o[0], dimb1, dimb0);
    out->element(k, a[0]) += ddot_(cden->size(), cden->data(), 1, cv0->data(), 1);
    out->element(k, a[1]) += ddot_(cden->size(), cden->data(), 1, cv1->data(), 1);
  }
  return out;
}

shared_ptr<GradFile> DKHcoreInfo::compute_v2(const array<shared_ptr<const Shell>,2>& s, const array<int,4>& a, const array<int,4>& o, const shared_ptr<const Matrix> den) const {
  assert(dkh2());
  const int dimb1 = s[0]->nbasis();
  const int dimb0 = s[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(o[1], o[0], dimb1, dimb0);
  
  auto out = make_shared<GradFile>(natom_);

  for (int k = 0; k < 3; k++) {
    std::shared_ptr<const Matrix> cv20 = v2relgrad_[3 * a[0] + k].get_submatrix(o[1], o[0], dimb1, dimb0);
    std::shared_ptr<const Matrix> cv21 = v2relgrad_[3 * a[1] + k].get_submatrix(o[1], o[0], dimb1, dimb0);
    out->element(k, a[0]) += ddot_(cden->size(), cden->data(), 1, cv20->data(), 1);
    out->element(k, a[1]) += ddot_(cden->size(), cden->data(), 1, cv21->data(), 1);
  }
  return out;
}


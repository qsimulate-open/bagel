//
// BAGEL - Parallel electron correlation program.
// Filename: moint.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/smith/moint.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


K2ext::K2ext(shared_ptr<const SMITH_Info> r, shared_ptr<const Coeff> c, vector<IndexRange> b) : ref_(r), coeff_(c), blocks_(b) {
  // so far MOInt can be called for 2-external K integral and all-internals.
  if (blocks_[0] != blocks_[2] || blocks_[1] != blocks_[3])
    throw logic_error("MOInt called with wrong blocks");
  data_ = make_shared<Tensor>(blocks_);
  form_4index(generate_list());
}


map<size_t, shared_ptr<DFFullDist>> K2ext::generate_list() {
  shared_ptr<const DFDist> df = ref_->geom()->df();

  // It is the easiest to do integral transformation for each blocks.
  assert(blocks_.size() == 4);
  map<size_t, shared_ptr<DFFullDist>> dflist;
  // AO dimension
  assert(df->nbasis0() == df->nbasis1());

  // occ loop
  for (auto& i0 : blocks_[0]) {
    shared_ptr<DFHalfDist> df_half = df->compute_half_transform(coeff_->slice(i0.offset(), i0.offset()+i0.size()))->apply_J();
    // virtual loop
    for (auto& i1 : blocks_[1]) {
      shared_ptr<DFFullDist> df_full = df_half->compute_second_transform(coeff_->slice(i1.offset(), i1.offset()+i1.size()));
      dflist.emplace(generate_hash_key(i0, i1), df_full);
    }
  }
  return dflist;
}


void K2ext::form_4index(const map<size_t, shared_ptr<DFFullDist>>& dflist) {
  // form four-index integrals
  // TODO this part should be heavily parallelized
  for (auto& i0 : blocks_[0]) {
    for (auto& i1 : blocks_[1]) {
      // find three-index integrals
      auto iter01 = dflist.find(generate_hash_key(i0, i1));
      assert(iter01 != dflist.end());
      shared_ptr<DFFullDist> df01 = iter01->second;
      size_t hashkey01 = generate_hash_key(i0, i1);

      for (auto& i2 : blocks_[2]) {
        for (auto& i3 : blocks_[3]) {
          // find three-index integrals
          size_t hashkey23 = generate_hash_key(i2, i3);
          if (hashkey23 > hashkey01) continue;

          auto iter23 = dflist.find(generate_hash_key(i2, i3));
          assert(iter23 != dflist.end());
          shared_ptr<const DFFullDist> df23 = iter23->second;

          // contract
          // TODO form_4index function now generates global 4 index tensor. This should be localized.
          shared_ptr<Matrix> tmp = df01->form_4index(df23, 1.0);
          unique_ptr<double[]> target(new double[tmp->size()]);
          copy_n(tmp->data(), tmp->size(), target.get()); // unnecessary copy

          // move in place
          if (hashkey23 != hashkey01) {
            unique_ptr<double[]> target2(new double[i0.size()*i1.size()*i2.size()*i3.size()]);
            blas::transpose(target.get(), i0.size()*i1.size(), i2.size()*i3.size(), target2.get());

            data_->put_block(target2, i2, i3, i0, i1);
          }

          data_->put_block(target, i0, i1, i2, i3);
        }
      }
    }
  }
}


MOFock::MOFock(shared_ptr<const SMITH_Info> r, vector<IndexRange> b) : ref_(r), coeff_(make_shared<Coeff>(*ref_->coeff())), blocks_(b) {
  // for simplicity, I assume that the Fock matrix is formed at once (may not be needed).
  assert(b.size() == 2 && b[0] == b[1]);
  const int ncore   = ref_->ncore();
  const int nclosed = ref_->nclosed() - ncore;
  assert(nclosed >= 0);
  const int nocc    = ref_->nocc();
  const int nact    = ref_->nact();
  const int nvirt   = ref_->nvirt();
  const int nbasis  = coeff_->ndim();

  data_  = make_shared<Tensor>(blocks_);
  h1_    = make_shared<Tensor>(blocks_);

  // cfock
  shared_ptr<Matrix> cfock = ref_->hcore()->copy();
  core_energy_ = 0.0;
  if (ncore+nclosed) {
    cfock = make_shared<Fock<1>>(r->geom(), ref_->hcore(), nullptr, coeff_->slice(0, ncore+nclosed), false, true);
    shared_ptr<const Matrix> den = coeff_->form_density_rhf(ncore+nclosed);
    core_energy_ = (*den * (*ref_->hcore()+*cfock)).trace() * 0.5;
  }

  shared_ptr<const Matrix> fock1;
  if (nact) {
    Matrix tmp(nact, nact);
    copy_n(ref_->rdm1_av()->data(), tmp.size(), tmp.data());
    tmp.sqrt();
    tmp.scale(1.0/sqrt(2.0));
    shared_ptr<Matrix> weighted_coeff = coeff_->slice_copy(ncore+nclosed, nocc);
    *weighted_coeff *= tmp;
    fock1 = make_shared<Fock<1>>(r->geom(), cfock, nullptr, weighted_coeff, false, true);
  } else {
    fock1 = cfock;
  }

  // We substitute diagonal part of the two-body integrals.
  // Note that E_ij,kl = E_ij E_kl - delta_jk E_il
  // and SMITH uses E_ij E_kl type excitations throughout. j and k must be active
  if (nact) {
    shared_ptr<const DFHalfDist> half = r->geom()->df()->compute_half_transform(coeff_->slice(ncore+nclosed, nocc))->apply_J();
    *cfock -= *half->form_2index(half, 0.5);
  }

  // if closed/virtual orbitals are present, we diagonalize the fock operator within this subspace
  const Matrix forig = *coeff_ % *fock1 * *coeff_;
  VectorB eig(nbasis);
  if (nclosed > 1) {
    shared_ptr<Matrix> fcl = forig.get_submatrix(ncore, ncore, nclosed, nclosed);
    fcl->diagonalize(eig);
    coeff_->copy_block(0, ncore, nbasis, nclosed, coeff_->slice(ncore, ncore+nclosed) * *fcl);
  }
  if (nvirt > 1) {
    shared_ptr<Matrix> fvirt = forig.get_submatrix(nocc, nocc, nvirt, nvirt);
    fvirt->diagonalize(eig);
    coeff_->copy_block(0, nocc, nbasis, nvirt, coeff_->slice(nocc, nocc+nvirt) * *fvirt);
  }
  const Matrix f = *coeff_ % *fock1 * *coeff_;
  const Matrix h1 = *coeff_ % *cfock * *coeff_;

  for (auto& i0 : blocks_[0]) {
    for (auto& i1 : blocks_[1]) {
      {
        shared_ptr<const Matrix> target = f.get_submatrix(i1.offset(), i0.offset(), i1.size(), i0.size());
        unique_ptr<double[]> tmp(new double[target->size()]);
        copy_n(target->data(), target->size(), tmp.get());
        data_->put_block(tmp, i1, i0);
      } {
        shared_ptr<const Matrix> target = h1.get_submatrix(i1.offset(), i0.offset(), i1.size(), i0.size());
        unique_ptr<double[]> tmp(new double[target->size()]);
        copy_n(target->data(), target->size(), tmp.get());
        h1_->put_block(tmp, i1, i0);
      }
    }
  }
}

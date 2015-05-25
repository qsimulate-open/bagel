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
#include <src/smith/smith_util.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/df/reldffull.h>
#include <src/scf/dhf/dfock.h>
#include <src/util/math/quatmatrix.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


template<typename DataType>
K2ext<DataType>::K2ext(shared_ptr<const SMITH_Info<DataType>> r, shared_ptr<const MatType> c, const vector<IndexRange>& b)
  : ref_(r), coeff_(c), blocks_(b) {

  // so far MOInt can be called for 2-external K integral and all-internals.
  if (blocks_[0] != blocks_[2] || blocks_[1] != blocks_[3])
    throw logic_error("MOInt called with wrong blocks");
  data_ = make_shared<Tensor_<DataType>>(blocks_);
  init();
}


template<>
void K2ext<complex<double>>::init() {
  // (1) make DFDists
  vector<shared_ptr<const DFDist>> dfs = ref_->geom()->dfs()->split_blocks();
  dfs.push_back(ref_->geom()->df());
  list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

  map<size_t, shared_ptr<RelDFFull>> dflist;

  // (2) first-transform
  for (auto& i0 : blocks_[0]) {
    list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, coeff_->slice_copy(i0.offset(), i0.offset()+i0.size()));
    for (auto& i : half_complex)
      i = i->apply_J();

    // (3) split and factorize
    list<shared_ptr<RelDFHalf>> half_complex_exch;
    for (auto& i : half_complex) {
      list<shared_ptr<RelDFHalf>> tmp = i->split(false);
      half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
    }
    half_complex.clear();
    DFock::factorize(half_complex_exch);

    for (auto& i1 : blocks_[1]) {
      // (4) compute (gamma|ia)
      list<shared_ptr<RelDFFull>> dffull;
      for (auto& i : half_complex_exch)
        dffull.push_back(make_shared<RelDFFull>(i, coeff_->slice_copy(i1.offset(), i1.offset()+i1.size())));
      DFock::factorize(dffull);
      dffull.front()->scale(dffull.front()->fac()); // take care of the factor
      assert(dffull.size() == 1);
      // adding this to dflist
      dflist.emplace(generate_hash_key(i0, i1), dffull.front());
    }
  }

  // form four-index integrals
  // TODO this part should be heavily parallelized
  for (auto& i0 : blocks_[0]) {
    for (auto& i1 : blocks_[1]) {
      // find three-index integrals
      auto iter01 = dflist.find(generate_hash_key(i0, i1));
      assert(iter01 != dflist.end());
      shared_ptr<RelDFFull> df01 = iter01->second;
      size_t hashkey01 = generate_hash_key(i0, i1);

      for (auto& i2 : blocks_[2]) {
        for (auto& i3 : blocks_[3]) {
          // find three-index integrals
          size_t hashkey23 = generate_hash_key(i2, i3);
          if (hashkey23 > hashkey01) continue;

          auto iter23 = dflist.find(generate_hash_key(i2, i3));
          assert(iter23 != dflist.end());
          shared_ptr<const RelDFFull> df23 = iter23->second;

          // contract
          // TODO form_4index function now generates global 4 index tensor. This should be localized.
          // conjugating because (ai|ai) is associated with an excitation operator
          shared_ptr<ZMatrix> tmp = df01->form_4index(df23, 1.0)->get_conjg();
          unique_ptr<complex<double>[]> target(new complex<double>[tmp->size()]);
          copy_n(tmp->data(), tmp->size(), target.get()); // unnecessary copy

          // move in place
          if (hashkey23 != hashkey01) {
            unique_ptr<complex<double>[]> target2(new complex<double>[i0.size()*i1.size()*i2.size()*i3.size()]);
            blas::transpose(target.get(), i0.size()*i1.size(), i2.size()*i3.size(), target2.get());
            data_->put_block(target2, i2, i3, i0, i1);
          }
          data_->put_block(target, i0, i1, i2, i3);
        }
      }
    }
  }
}


template<>
void K2ext<double>::init() {
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


template<typename DataType>
MOFock<DataType>::MOFock(shared_ptr<const SMITH_Info<DataType>> r, const vector<IndexRange>& b) : ref_(r), coeff_(ref_->coeff()), blocks_(b) {
  assert(b.size() == 2 && b[0] == b[1]);

  data_  = make_shared<Tensor_<DataType>>(blocks_);
  h1_    = make_shared<Tensor_<DataType>>(blocks_);
  init();
}


template<>
void MOFock<complex<double>>::init() {
  const int ncore   = ref_->ncore();
  const int nclosed = ref_->nclosed() - ncore;
  assert(nclosed >= 0);
  const int nocc    = ref_->nocc();
  const int nact    = ref_->nact();
  const int nvirt   = ref_->nvirt();

  auto relref = dynamic_pointer_cast<const RelReference>(ref_->ref());

  // first hcore
  auto hcore = make_shared<RelHcore>(ref_->geom());
  shared_ptr<ZMatrix> cfock;
  core_energy_ = 0.0;
  if (ncore+nclosed) {
    cfock = make_shared<DFock>(ref_->geom(), hcore, coeff_->slice_copy(0, 2*(ncore+nclosed)), /*gaunt*/false, /*breit*/false, /*store_half*/false);
    shared_ptr<const ZMatrix> den = coeff_->form_density_rhf(2*(ncore+nclosed), 0);
    core_energy_ = detail::real((*den * (*hcore+*cfock)).trace()) * 0.5;
  } else {
    cfock = hcore->copy();
  }

  // active fock
  shared_ptr<const ZMatrix> fock1;
  if (nact) {
    shared_ptr<ZMatrix> tmp = relref->rdm1_av()->get_conjg();
    tmp->sqrt();
    shared_ptr<ZMatrix> weighted_coeff = coeff_->slice_copy(2*(ncore+nclosed), 2*nocc);
    *weighted_coeff *= *tmp;
    fock1 = make_shared<DFock>(ref_->geom(), cfock, weighted_coeff, false, false, false);
  } else {
    fock1 = cfock;
  }

  // We substitute diagonal part of the two-body integrals.
  // Note that E_ij,kl = E_ij E_kl - delta_jk E_il
  // and SMITH uses E_ij E_kl type excitations throughout. j and k must be active
  if (nact) {
    vector<shared_ptr<const DFDist>> dfs = ref_->geom()->dfs()->split_blocks();
    dfs.push_back(ref_->geom()->df());
    list<shared_ptr<RelDF>> dfdists = DFock::make_dfdists(dfs, false);

    list<shared_ptr<RelDFHalf>> half_complex = DFock::make_half_complex(dfdists, coeff_->slice_copy(2*(ncore+nclosed), 2*nocc));
    for (auto& i : half_complex)
      i = i->apply_J();

    list<shared_ptr<RelDFHalf>> half_complex_exch;
    for (auto& i : half_complex) {
      list<shared_ptr<RelDFHalf>> tmp = i->split(false);
      half_complex_exch.insert(half_complex_exch.end(), tmp.begin(), tmp.end());
    }
    half_complex.clear();
    DFock::factorize(half_complex_exch);

    for (auto& i : half_complex_exch)
      i->set_sum_diff();

    // computing K operators
    int icnt = 0;
    for (auto& i : half_complex_exch) {
      int jcnt = 0;
      for (auto& j : half_complex_exch) {
        if (i->alpha_matches(j) && icnt <= jcnt)
          DFock::add_Exop_block(*cfock, i, j, 0.5, icnt == jcnt);
        ++jcnt;
      }
      ++icnt;
    }
  }

  // if closed/virtual orbitals are present, we diagonalize the fock operator within this subspace
  const ZMatrix forig = *coeff_ % *fock1 * *coeff_;
  VectorB eig(forig.ndim());
  auto newcoeff = coeff_->copy();
  if (nclosed > 1) {
    auto fcl = make_shared<QuatMatrix>(*forig.get_submatrix(0, 0, (ncore+nclosed)*2, (ncore+nclosed)*2));
    fcl->diagonalize(eig);
    newcoeff->copy_block(0, 0, newcoeff->ndim(), (ncore+nclosed)*2, newcoeff->slice(0, (ncore+nclosed)*2) * *fcl);
  }
  if (nvirt > 1) {
    auto fvirt = make_shared<QuatMatrix>(*forig.get_submatrix(nocc*2, nocc*2, nvirt*2, nvirt*2));
    fvirt->diagonalize(eig);
    newcoeff->copy_block(0, nocc*2, newcoeff->ndim(), nvirt*2, newcoeff->slice(nocc*2, (nocc+nvirt)*2) * *fvirt);
  }
  // **** CAUTION **** updating the coefficient
  coeff_ = newcoeff;

  auto f  = make_shared<ZMatrix>(*coeff_ % *fock1 * *coeff_);
  auto h1 = make_shared<ZMatrix>(*coeff_ % *cfock * *coeff_);

  if (!f->is_hermitian()) throw logic_error("Fock is not Hermitian");
  if (!h1->is_hermitian()) throw logic_error("Hcore is not Hermitian");

  fill_block<2,complex<double>>(data_, f->get_conjg(), {0,0}, blocks_);
  fill_block<2,complex<double>>(h1_,  h1->get_conjg(), {0,0}, blocks_);
}

template<>
void MOFock<double>::init() {
  // for simplicity, I assume that the Fock matrix is formed at once (may not be needed).
  const int ncore   = ref_->ncore();
  const int nclosed = ref_->nclosed() - ncore;
  assert(nclosed >= 0);
  const int nocc    = ref_->nocc();
  const int nact    = ref_->nact();
  const int nvirt   = ref_->nvirt();
  const int nbasis  = coeff_->ndim();

  // cfock
  shared_ptr<Matrix> cfock = ref_->hcore()->copy();
  core_energy_ = 0.0;
  if (ncore+nclosed) {
    cfock = make_shared<Fock<1>>(ref_->geom(), ref_->hcore(), nullptr, coeff_->slice(0, ncore+nclosed), false, true);
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
    fock1 = make_shared<Fock<1>>(ref_->geom(), cfock, nullptr, weighted_coeff, false, true);
  } else {
    fock1 = cfock;
  }

  // We substitute diagonal part of the two-body integrals.
  // Note that E_ij,kl = E_ij E_kl - delta_jk E_il
  // and SMITH uses E_ij E_kl type excitations throughout. j and k must be active
  if (nact) {
    shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_->slice(ncore+nclosed, nocc))->apply_J();
    *cfock -= *half->form_2index(half, 0.5);
  }

  // if closed/virtual orbitals are present, we diagonalize the fock operator within this subspace
  const Matrix forig = *coeff_ % *fock1 * *coeff_;
  VectorB eig(nbasis);
  auto newcoeff = coeff_->copy();
  if (nclosed > 1) {
    shared_ptr<Matrix> fcl = forig.get_submatrix(0, 0, ncore+nclosed, ncore+nclosed);
    fcl->diagonalize(eig);
    newcoeff->copy_block(0, 0, nbasis, ncore+nclosed, newcoeff->slice(0, ncore+nclosed) * *fcl);
  }
  if (nvirt > 1) {
    shared_ptr<Matrix> fvirt = forig.get_submatrix(nocc, nocc, nvirt, nvirt);
    fvirt->diagonalize(eig);
    newcoeff->copy_block(0, nocc, nbasis, nvirt, newcoeff->slice(nocc, nocc+nvirt) * *fvirt);
  }
  // **** CAUTION **** updating the coefficient
  coeff_ = newcoeff;
  auto f  = make_shared<Matrix>(*coeff_ % *fock1 * *coeff_);
  auto h1 = make_shared<Matrix>(*coeff_ % *cfock * *coeff_);

  fill_block<2,double>(data_,  f, {0,0}, blocks_);
  fill_block<2,double>(h1_,   h1, {0,0}, blocks_);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class K2ext<double>;
template class K2ext<complex<double>>;
template class MOFock<double>;
template class MOFock<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

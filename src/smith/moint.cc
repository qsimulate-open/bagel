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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/moint.h>
#include <src/smith/smith_util.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/df/reldffull.h>
#include <src/scf/dhf/dfock.h>
#include <src/ci/zfci/relmofile.h>
#include <src/util/math/quatmatrix.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


template<typename DataType>
K2ext<DataType>::K2ext(shared_ptr<const SMITH_Info<DataType>> r, shared_ptr<const MatType> c, const vector<IndexRange>& b)
  : info_(r), coeff_(c), blocks_(b) {

  // so far MOInt can be called for 2-external K integral and all-internals.
  if (blocks_[0] != blocks_[2] || blocks_[1] != blocks_[3])
    throw logic_error("MOInt called with wrong blocks");
  data_ = make_shared<TATensor<DataType,4>>(blocks_, false, is_same<DataType,complex<double>>::value);
  init();
}


template<>
void K2ext<complex<double>>::init() {

// broken ... waiting to hear symmetry support in TiledArray
#if 0
  // bits to store
  const bool braket = blocks_[0] == blocks_[1] && blocks_[2] == blocks_[3];
  const vector<vector<int>> cblocks = braket ?
    vector<vector<int>>{{0,0,0,0}, {0,0,0,1}, {0,0,1,1}, {0,1,0,1}, {0,1,1,0}, {0,1,1,1}, {1,1,1,1}} :
    vector<vector<int>>{{0,0,0,0}, {0,0,0,1}, {0,0,1,0}, {0,0,1,1}, {0,1,0,1}, {0,1,1,0}, {0,1,1,1}, {1,0,1,1}, {1,0,1,0}, {1,1,1,1}};

  auto compute = [this, &cblocks](const bool gaunt, const bool breit) {

    map<size_t, shared_ptr<const ListRelDFFull>> dflist, dflist2;

    for (auto& i0 : blocks_[0]) {
      shared_ptr<const ZMatrix> i0coeff = coeff_->slice_copy(i0.offset(), i0.offset()+i0.size());
      list<shared_ptr<RelDFHalf>> half, half2;
      tie(half, half2) = RelMOFile::compute_half(info_->geom(), i0coeff, gaunt, breit);

      for (auto& i1 : blocks_[1]) {
        shared_ptr<const ZMatrix> i1coeff = coeff_->slice_copy(i1.offset(), i1.offset()+i1.size());
        dflist.emplace(generate_hash_key(i0, i1), RelMOFile::compute_full(i1coeff, half, true));
        if (breit)
          dflist2.emplace(generate_hash_key(i0, i1), RelMOFile::compute_full(i1coeff, half2, false));
      }
    }
    if (!breit)
      dflist2 = dflist;

    // form four-index integrals
    // TODO this part should be heavily parallelized
    const double gscale = gaunt ? (breit ? -0.25 /*we explicitly symmetrize*/ : -1.0) : 1.0;
    for (auto& i0 : blocks_[0]) {
      for (auto& i1 : blocks_[1]) {
        // find three-index integrals
        size_t hashkey01 = generate_hash_key(i0, i1);
        assert(dflist.find(hashkey01) != dflist.end());
        shared_ptr<const ListRelDFFull> df01   = dflist.find(hashkey01)->second;
        shared_ptr<const ListRelDFFull> df01_2 = dflist2.find(hashkey01)->second;
        const int t0 = i0.kramers() ? 1 : 0;
        const int t1 = i1.kramers() ? 1 : 0;

        for (auto& i2 : blocks_[2]) {
          for (auto& i3 : blocks_[3]) {
            // find three-index integrals
            const int t2 = i2.kramers() ? 1 : 0;
            const int t3 = i3.kramers() ? 1 : 0;
            if (find(cblocks.begin(), cblocks.end(), vector<int>{t0, t1, t2, t3}) == cblocks.end())
              continue;

            size_t hashkey23 = generate_hash_key(i2, i3);
            assert(dflist.find(hashkey23) != dflist.end());
            shared_ptr<const ListRelDFFull> df23   = dflist.find(hashkey23)->second;
            shared_ptr<const ListRelDFFull> df23_2 = dflist2.find(hashkey23)->second;

            // contract
            // TODO form_4index function now generates global 4 index tensor. This should be localized.
            // conjugating because (ai|ai) is associated with an excitation operator
            unique_ptr<complex<double>[]> target = data_->move_block(i0, i1, i2, i3);
            {
              shared_ptr<ZMatrix> tmp = df01->form_4index(df23_2, 1.0)->get_conjg();
              blas::ax_plus_y_n(gscale, tmp->data(), tmp->size(), target.get());
            }
            if (breit) {
              shared_ptr<ZMatrix> tmp = df01_2->form_4index(df23, 1.0)->get_conjg();
              blas::ax_plus_y_n(gscale, tmp->data(), tmp->size(), target.get());
            }
            data_->put_block(target, i0, i1, i2, i3);
          }
        }
      }
    }
  };

  // coulomb operator
  compute(false, false);
  if (info_->gaunt())
    compute(true, info_->breit());

  map<vector<int>, pair<double,bool>> perm{{{0,1,2,3}, {1.0, false}}, {{2,3,0,1}, {1.0, false}}};
  if (braket) {
    perm.emplace(vector<int>{1,0,3,2}, make_pair(1.0, true));
    perm.emplace(vector<int>{3,2,1,0}, make_pair(1.0, true));
  }
  data_->set_perm(perm);
#endif
}


template<>
void K2ext<double>::init() {
  shared_ptr<const DFDist> df = info_->geom()->df();

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
          shared_ptr<Matrix> tmp = df01->form_4index(df23, 1.0);

          {
            const std::vector<Index> index = {i0, i1, i2, i3};
            auto tileiter = data_->get_local(index);
            if (tileiter.first) {
              const TiledArray::Range range = data_->trange().make_tile_range(tileiter.second.ordinal());
              typename TiledArray::Array<double,4>::value_type tile(range);
              copy_n(tmp->data(), tmp->size(), &(tile[0]));
              *tileiter.second = tile;
            }
          }

          if (hashkey23 != hashkey01) {
            const std::vector<Index> index = {i2, i3, i0, i1};
            auto tileiter = data_->get_local(index);
            if (tileiter.first) {
              const TiledArray::Range range = data_->trange().make_tile_range(tileiter.second.ordinal());
              typename TiledArray::Array<double,4>::value_type tile(range);
              blas::transpose(tmp->data(), i0.size()*i1.size(), i2.size()*i3.size(), &(tile[0]));
              *tileiter.second = tile;
            }
          }
        }
      }
    }
  }
}


template<typename DataType>
MOFock<DataType>::MOFock(shared_ptr<const SMITH_Info<DataType>> r, const vector<IndexRange>& b) : info_(r), coeff_(info_->coeff()), blocks_(b) {
  assert(b.size() == 2 && b[0] == b[1]);

  data_  = make_shared<TATensor<DataType,2>>(blocks_);
  h1_    = make_shared<TATensor<DataType,2>>(blocks_);
  init();
}


template<>
void MOFock<complex<double>>::init() {
  const int ncore   = info_->ncore();
  const int nclosed = info_->nclosed() - ncore;
  assert(nclosed >= 0);
  const int nocc    = info_->nocc();
  const int nact    = info_->nact();

  auto relref = dynamic_pointer_cast<const RelReference>(info_->ref());

  // first hcore
  auto hcore = make_shared<RelHcore>(info_->geom());
  shared_ptr<ZMatrix> cfock;
  core_energy_ = 0.0;
  if (ncore+nclosed) {
    cfock = make_shared<DFock>(info_->geom(), hcore, coeff_->slice_copy(0, 2*(ncore+nclosed)),
                               info_->gaunt(), info_->breit(), /*store_half*/false, info_->breit());
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
    fock1 = make_shared<DFock>(info_->geom(), cfock, weighted_coeff, info_->gaunt(), info_->breit(), false, info_->breit());
  } else {
    fock1 = cfock;
  }

  // We substitute diagonal part of the two-body integrals.
  // Note that E_ij,kl = E_ij E_kl - delta_jk E_il
  // and SMITH uses E_ij E_kl type excitations throughout. j and k must be active
  if (nact) {
    cfock = make_shared<DFock>(info_->geom(), cfock, coeff_->slice_copy(2*(ncore+nclosed), 2*nocc),
                               info_->gaunt(), info_->breit(), /*store_half*/false, /*robust*/info_->breit(),
                               /*scale_exch*/0.5, /*scale_coulomb*/0.0);
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
  const int nvirt = info_->nvirt();
  const int nvirtall = nvirt+info_->nfrozenvirt();
  if (nvirtall > 1) {
    auto fvirt = make_shared<QuatMatrix>(*forig.get_submatrix(nocc*2, nocc*2, nvirtall*2, nvirtall*2));
    fvirt->diagonalize(eig);
    const ZMatrix crot = coeff_->slice(nocc*2, (nocc+nvirtall)*2) * *fvirt;
    newcoeff->copy_block(0, nocc*2,       newcoeff->ndim(), nvirt, crot.slice(0, nvirt));
    newcoeff->copy_block(0, nocc*2+nvirt, newcoeff->ndim(), nvirt, crot.slice(nvirtall, nvirtall+nvirt));
    if (info_->nfrozenvirt() > 0) {
      cout << "       - Truncating virtual orbitals: " << setw(20) << setprecision(10) << eig[nvirt] << endl;
      newcoeff = newcoeff->slice_copy(0, (nocc+nvirt)*2);
    }
  }

  // **** CAUTION **** updating the coefficient
  coeff_ = newcoeff;

  auto f  = make_shared<ZMatrix>(*coeff_ % *fock1 * *coeff_);
  auto h1 = make_shared<ZMatrix>(*coeff_ % *cfock * *coeff_);

  if (!f->is_hermitian()) throw logic_error("Fock is not Hermitian");
  if (!h1->is_hermitian()) throw logic_error("Hcore is not Hermitian");

  {
    unique_ptr<complex<double>[]> eig = f->diag();
    const int n = f->ndim();
    eig_ = vector<double>(n);
    for (int i = 0; i != n; ++i) {
      assert(fabs(imag(eig[i])) < 1.0e-10);
      eig_[i] = real(eig[i]);
    }
  }

  fill_block<2,complex<double>>(data_, f->get_conjg(), {0,0});
  fill_block<2,complex<double>>(h1_,  h1->get_conjg(), {0,0});
}

template<>
void MOFock<double>::init() {
  // for simplicity, I assume that the Fock matrix is formed at once (may not be needed).
  const int ncore   = info_->ncore();
  const int nclosed = info_->nclosed() - ncore;
  assert(nclosed >= 0);
  const int nocc    = info_->nocc();
  const int nact    = info_->nact();
  const int nbasis  = coeff_->ndim();

  // cfock
  shared_ptr<Matrix> cfock = info_->hcore()->copy();
  core_energy_ = 0.0;
  if (ncore+nclosed) {
    cfock = make_shared<Fock<1>>(info_->geom(), info_->hcore(), nullptr, coeff_->slice(0, ncore+nclosed), false, true);
    shared_ptr<const Matrix> den = coeff_->form_density_rhf(ncore+nclosed);
    core_energy_ = (*den * (*info_->hcore()+*cfock)).trace() * 0.5;
  }

  shared_ptr<const Matrix> fock1;
  if (nact) {
    Matrix tmp(nact, nact);
    copy_n(info_->rdm1_av()->data(), tmp.size(), tmp.data());
    tmp.sqrt();
    tmp.scale(1.0/sqrt(2.0));
    shared_ptr<Matrix> weighted_coeff = coeff_->slice_copy(ncore+nclosed, nocc);
    *weighted_coeff *= tmp;
    fock1 = make_shared<Fock<1>>(info_->geom(), cfock, nullptr, weighted_coeff, false, true);
  } else {
    fock1 = cfock;
  }

  // We substitute diagonal part of the two-body integrals.
  // Note that E_ij,kl = E_ij E_kl - delta_jk E_il
  // and SMITH uses E_ij E_kl type excitations throughout. j and k must be active
  if (nact) {
    shared_ptr<const DFHalfDist> half = info_->geom()->df()->compute_half_transform(coeff_->slice(ncore+nclosed, nocc))->apply_J();
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
  const int nvirt   = info_->nvirt();
  const int nvirtall = nvirt+info_->nfrozenvirt();
  if (nvirtall > 1) {
    shared_ptr<Matrix> fvirt = forig.get_submatrix(nocc, nocc, nvirtall, nvirtall);
    fvirt->diagonalize(eig);
    newcoeff->copy_block(0, nocc, nbasis, nvirtall, newcoeff->slice(nocc, nocc+nvirtall) * *fvirt);
    if (info_->nfrozenvirt() > 0) {
      cout << "       - Truncating virtual orbitals: " << setw(20) << setprecision(10) << eig[nvirt] << endl;
      newcoeff = newcoeff->slice_copy(0, nocc+nvirt);
    }
  }

  // **** CAUTION **** updating the coefficient
  coeff_ = newcoeff;
  auto f  = make_shared<Matrix>(*coeff_ % *fock1 * *coeff_);
  auto h1 = make_shared<Matrix>(*coeff_ % *cfock * *coeff_);

  {
    unique_ptr<double[]> eig = f->diag();
    const int n = f->ndim();
    eig_ = vector<double>(eig.get(), eig.get()+n);
  }

  fill_block<2,double>(data_,  f, {0,0});
  fill_block<2,double>(h1_,   h1, {0,0});
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class K2ext<double>;
template class K2ext<complex<double>>;
template class MOFock<double>;
template class MOFock<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

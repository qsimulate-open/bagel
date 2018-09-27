//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: moint.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/moint.h>
#include <src/smith/smith_util.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/df/reldffull.h>
#include <src/scf/dhf/dfock.h>
#include <src/ci/zfci/reljop.h>
#include <src/util/math/quatmatrix.h>
#include <src/util/exception.h>
#include <iostream>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


template<typename DataType>
K2ext<DataType>::K2ext(shared_ptr<const SMITH_Info<DataType>> r, shared_ptr<const MatType> c, const vector<IndexRange>& b)
  : info_(r), coeff_(c), blocks_(b) {

  // so far MOInt can be called for 2-external K integral and all-internals.
  if (blocks_[0] != blocks_[2] || blocks_[1] != blocks_[3])
    throw logic_error("MOInt called with wrong blocks");
  init();
}


template<>
void K2ext<complex<double>>::init() {
  // bits to store
  const bool braket = blocks_[0] == blocks_[1] && blocks_[2] == blocks_[3];
  const list<vector<int>> cblocks_int = braket ?
    list<vector<int>>{{0,0,0,0}, {0,0,0,1}, {0,0,1,1}, {0,1,0,1}, {0,1,1,0}, {0,1,1,1}, {1,1,1,1}} :
    list<vector<int>>{{0,0,0,0}, {0,0,0,1}, {0,0,1,0}, {0,0,1,1}, {0,1,0,1}, {0,1,1,0}, {0,1,1,1}, {1,0,1,1}, {1,0,1,0}, {1,1,1,1}};
  // convert this to list<vector<bool>>
  list<vector<bool>> cblocks;
  for (auto& i : cblocks_int) {
    vector<bool> tmp;
    for (auto& j : i)
      tmp.push_back(static_cast<bool>(j));
    cblocks.push_back(tmp);
  }
  // tabluate all the blocks that are to be stored
  unordered_set<size_t> sparse;
  for (auto& i0 : blocks_[0])
    for (auto& i1 : blocks_[1])
      for (auto& i2 : blocks_[2])
        for (auto& i3 : blocks_[3])
          if (find(cblocks.begin(), cblocks.end(), vector<bool>{i0.kramers(), i1.kramers(), i2.kramers(), i3.kramers()}) != cblocks.end())
            sparse.insert(generate_hash_key(i0, i1, i2, i3));
  // allocate
  data_ = make_shared<Tensor_<complex<double>>>(blocks_, /*kramers*/true, sparse, /*alloc*/true);
  data_->set_stored_sectors(cblocks);

  // Aux index blocking
  const IndexRange aux(info_->geom()->df()->adist_now());
  const size_t astart = info_->geom()->df()->block(0)->astart();

  auto compute = [&, this](const bool gaunt, const bool breit) {
    // create an intermediate array
    using MapType = map<int, shared_ptr<Tensor_<complex<double>>>>;
    MapType ext, ext2;
    auto alpha = gaunt ? list<int>{Comp::X, Comp::Y, Comp::Z} : list<int>{Comp::L};

    for (auto& i : alpha) {
      ext.emplace(i, make_shared<Tensor_<complex<double>>>(vector<IndexRange>{aux, blocks_[0], blocks_[1]}, false, unordered_set<size_t>{}, true));
      if (breit)
        ext2.emplace(i, make_shared<Tensor_<complex<double>>>(vector<IndexRange>{aux, blocks_[0], blocks_[1]}, false, unordered_set<size_t>{}, true));
    }

    for (auto& i0 : blocks_[0]) {
      shared_ptr<const ZMatrix> i0coeff = coeff_->slice_copy(i0.offset(), i0.offset()+i0.size());
      list<shared_ptr<RelDFHalf>> half, half2;
      tie(half, half2) = RelJop::compute_half(info_->geom(), i0coeff, gaunt, breit);

      for (auto& i1 : blocks_[1]) {
        shared_ptr<const ZMatrix> i1coeff = coeff_->slice_copy(i1.offset(), i1.offset()+i1.size());
        shared_ptr<const ListRelDFFull> full, full2;
        full = RelJop::compute_full(i1coeff, half, true);
        if (breit)
          full2 = RelJop::compute_full(i1coeff, half2, false);

        for (auto& a : aux)
          if (a.offset() == astart) {
            const size_t bufsize = a.size()*i0.size()*i1.size();
            unique_ptr<complex<double>[]> buf(new complex<double>[bufsize]);
            auto comp = [&](shared_ptr<const ListRelDFFull> cfull, MapType& cext) {
              for (auto& data : cfull->data()) {
                auto block = data->get_block(a.offset(), a.size(), 0, i0.size(), 0, i1.size());
                assert(block->size() == bufsize);
                copy_n(block->data(), bufsize, buf.get());
                cext.at(data->alpha_comp())->put_block(buf, a, i0, i1);
              }
            };
            comp(full, ext);
            if (breit)
              comp(full2, ext2);
          }
      }
    }
    if (!breit) ext2 = ext;
    // wait for other nodes
    for (auto& i : ext)
      i.second->fence();
    for (auto& i : ext2)
      i.second->fence();

    // form four-index integrals
    const double gscale = gaunt ? (breit ? -0.25 /*we explicitly symmetrize*/ : -1.0) : 1.0;
    for (auto& i0 : blocks_[0]) {
      for (auto& i1 : blocks_[1]) {
        for (auto& i2 : blocks_[2]) {
          for (auto& i3 : blocks_[3]) {
            if (sparse.count(generate_hash_key(i0, i1, i2, i3)) == 0 || !data_->is_local(i0, i1, i2, i3))
              continue;
            const size_t bufsize = data_->get_size(i0, i1, i2, i3);
            unique_ptr<complex<double>[]> buf(new complex<double>[bufsize]);
            fill_n(buf.get(), bufsize, 0.0);

            for (auto& a : aux) {
              for (auto& i : alpha) {
                auto comp = [&] (const MapType& cext, const MapType& cext2) {
                  unique_ptr<complex<double>[]> data01 = cext.at(i)->get_block(a, i0, i1);
                  unique_ptr<complex<double>[]> data23 = cext2.at(i)->get_block(a, i2, i3);
                  btas::gemm_impl<true>::call(CblasColMajor, CblasTrans, CblasNoTrans, i0.size()*i1.size(), i2.size()*i3.size(), a.size(),
                                              gscale, data01.get(), a.size(), data23.get(), a.size(), 1.0, buf.get(), i0.size()*i1.size());
                };
                comp(ext, ext2);
                if (breit)
                  comp(ext2, ext);
              }
            }
            blas::conj_n(buf.get(), bufsize);
            data_->add_block(buf, i0, i1, i2, i3);
          }
        }
      }
    }
    data_->fence();
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
}


template<>
void K2ext<double>::init() {
  data_ = make_shared<Tensor_<double>>(blocks_);
  data_->allocate();

  shared_ptr<const DFDist> df = info_->geom()->df();

  // It is the easiest to do integral transformation for each blocks.
  assert(blocks_.size() == 4);
  // AO dimension
  assert(df->nbasis0() == df->nbasis1());

  // Aux index blocking
  const IndexRange aux(df->adist_now());

  // create an intermediate array
  Tensor_<double> ext(vector<IndexRange>{aux, blocks_[0], blocks_[1]});
  ext.allocate();

  // occ loop
  for (auto& i0 : blocks_[0]) {
    shared_ptr<DFHalfDist> df_half = df->compute_half_transform(coeff_->slice(i0.offset(), i0.offset()+i0.size()))->apply_J();
    // virtual loop
    for (auto& i1 : blocks_[1]) {
      shared_ptr<DFFullDist> df_full = df_half->compute_second_transform(coeff_->slice(i1.offset(), i1.offset()+i1.size()));
      const size_t bufsize = df_full->block(0)->size();
      unique_ptr<double[]> buf(new double[bufsize]);
      copy_n(df_full->block(0)->data(), bufsize, buf.get());

      for (auto& a : aux)
        if (a.offset() == df->block(0)->astart())
          ext.put_block(buf, a, i0, i1);
    }
  }
  // wait for other nodes
  ext.fence();

  // form four-index integrals
  for (auto& i0 : blocks_[0]) {
    for (auto& i1 : blocks_[1]) {
      for (auto& i2 : blocks_[2]) {
        for (auto& i3 : blocks_[3]) {
          const size_t hashkey01 = generate_hash_key(i0, i1);
          const size_t hashkey23 = generate_hash_key(i2, i3);
          if (hashkey23 > hashkey01) continue;
          if (!data_->is_local(i0, i1, i2, i3)) continue;

          const size_t bufsize = data_->get_size(i0, i1, i2, i3);
          unique_ptr<double[]> buf0(new double[bufsize]);
          fill_n(buf0.get(), bufsize, 0.0);

          for (auto& a : aux) {
            unique_ptr<double[]> data01 = ext.get_block(a, i0, i1);
            unique_ptr<double[]> data23 = ext.get_block(a, i2, i3);
            // contract and accumulate
            btas::gemm_impl<true>::call(CblasColMajor, CblasTrans, CblasNoTrans, i0.size()*i1.size(), i2.size()*i3.size(), a.size(),
                                        1.0, data01.get(), a.size(), data23.get(), a.size(), 1.0, buf0.get(), i0.size()*i1.size());
          }

          // put in place
          data_->put_block(buf0, i0, i1, i2, i3);

          if (hashkey23 != hashkey01) {
            unique_ptr<double[]> buf1(new double[bufsize]);
            blas::transpose(buf0.get(), i0.size()*i1.size(), i2.size()*i3.size(), buf1.get());
            data_->put_block(buf1, i2, i3, i0, i1);
          }
        }
      }
    }
  }
  data_->fence();
}


template<typename DataType>
MOFock<DataType>::MOFock(shared_ptr<const SMITH_Info<DataType>> r, const vector<IndexRange>& b) : info_(r), coeff_(info_->coeff()), blocks_(b) {
  assert(b.size() == 2 && b[0] == b[1]);

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
    assert(fcl->is_t_symmetric(1.0e-6));
    fcl->diagonalize(eig);
    newcoeff->copy_block(0, 0, newcoeff->ndim(), (ncore+nclosed)*2, newcoeff->slice(0, (ncore+nclosed)*2) * *fcl);
    if (ncore)
      cout << endl << "  Core orbital energies: " << endl;
    for (int i = 0; i != ncore+nclosed; ++i) {
      if (i == ncore)
        cout << endl << "  Closed orbital energies: " << endl;
      cout << "    " << i+1 << "  " << setw(12) << setprecision(4) << fixed << eig[i] << endl;
    }
    cout << endl;
  }

  const int nvirt = info_->nvirt();
  const int nvirtall = nvirt+info_->nfrozenvirt();
  if (nvirtall > 1) {
    auto fvirt = make_shared<QuatMatrix>(*forig.get_submatrix(nocc*2, nocc*2, nvirtall*2, nvirtall*2));
    assert(fvirt->is_t_symmetric(1.0e-6));
    fvirt->diagonalize(eig);

    const ZMatrix crot = coeff_->slice(nocc*2, (nocc+nvirtall)*2) * *fvirt;
    newcoeff->copy_block(0, nocc*2,       newcoeff->ndim(), nvirt, crot.slice(0, nvirt));
    newcoeff->copy_block(0, nocc*2+nvirt, newcoeff->ndim(), nvirt, crot.slice(nvirtall, nvirtall+nvirt));
    if (info_->nfrozenvirt() > 0) {
      cout << "       - Truncating virtual orbitals: " << setw(20) << setprecision(10) << eig[nvirt] << endl;
      newcoeff = newcoeff->slice_copy(0, (nocc+nvirt)*2);
    }
    cout << endl << "  Virtual orbital energies: " << endl;
    for (int i = 0; i != nvirtall; ++i) {
      if (i == nvirt)
        cout << endl << "  Deleted orbital energies: " << endl;
      cout << "    " << nvirtall - i << "  " << setw(12) << setprecision(4) << fixed << eig[i] << endl;
    }
    cout << endl;
  }

  // **** CAUTION **** updating the coefficient
  coeff_ = newcoeff;

  auto f  = make_shared<ZMatrix>(*coeff_ % *fock1 * *coeff_);
  auto h1 = make_shared<ZMatrix>(*coeff_ % *cfock * *coeff_);

  if (!f->is_hermitian()) throw logic_error("Fock is not Hermitian");
  if (!h1->is_hermitian()) throw logic_error("Hcore is not Hermitian");

  if (info_->block_diag_fock()) {
    cout << "  * Removing off-diagonal blocks of the relativistic Fock matrix" << endl;
    if (to_lower(info_->method()) == "casa")
      cout << "    CAS/A with these blocks neglected is equivalent to partially contracted NEVPT2." << endl;
    auto fsave = f->copy();
    f->zero();
    const int nc = 2;
    assert(f->ndim() == f->mdim() && f->ndim() == nc * info_->nocc() + nc * info_->nvirt());
    f->copy_block(0, 0, nc*info_->nclosed(), nc*info_->nclosed(), fsave->get_submatrix(0, 0, nc*info_->nclosed(), nc*info_->nclosed()));
    f->copy_block(nc*info_->nclosed(), nc*info_->nclosed(), nc*info_->nact(), nc*info_->nact(),
                  fsave->get_submatrix(nc*info_->nclosed(), nc*info_->nclosed(), nc*info_->nact(), nc*info_->nact()));
    f->copy_block(nc*info_->nocc(), nc*info_->nocc(), nc*info_->nvirt(), nc*info_->nvirt(),
                  fsave->get_submatrix(nc*info_->nocc(), nc*info_->nocc(), nc*info_->nvirt(), nc*info_->nvirt()));
  }

  data_ = fill_block<2,complex<double>>(f->get_conjg(), {0,0}, blocks_);
  h1_   = fill_block<2,complex<double>>(h1->get_conjg(), {0,0}, blocks_);

  // print out the Fock matrix and terminate when preparing for external RDM runs
  if (info_->external_rdm() == "noref") {
    stringstream ss; ss << scientific << setprecision(15);
    for (int i = 0; i != nact*2; ++i) {
      for (int j = 0; j != nact*2; ++j)
        ss << f->element(j+(nclosed+ncore)*2, i+(nclosed+ncore)*2) << " ";
      ss << endl;
    }
    ofstream fs("FOCKMAT");
    fs << ss.str();
    throw Termination("Fock matrix has been dumped in FOCKMAT");
  }
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

  if (info_->block_diag_fock()) {
    cout << "  * Removing off-diagonal blocks of the (nonrel) Fock matrix" << endl;
    if (to_lower(info_->method()) == "casa")
      cout << "    CAS/A with these blocks neglected is equivalent to partially contracted NEVPT2." << endl;
    auto fsave = f->copy();
    f->zero();
    const int nc = 1;
    assert(f->ndim() == f->mdim() && f->ndim() == nc * info_->nocc() + nc * info_->nvirt());
    f->copy_block(0, 0, nc*info_->nclosed(), nc*info_->nclosed(), fsave->get_submatrix(0, 0, nc*info_->nclosed(), nc*info_->nclosed()));
    f->copy_block(nc*info_->nclosed(), nc*info_->nclosed(), nc*info_->nact(), nc*info_->nact(), fsave->get_submatrix(nc*info_->nclosed(), nc*info_->nclosed(), nc*info_->nact(), nc*info_->nact()));
    f->copy_block(nc*info_->nocc(), nc*info_->nocc(), nc*info_->nvirt(), nc*info_->nvirt(), fsave->get_submatrix(nc*info_->nocc(), nc*info_->nocc(), nc*info_->nvirt(), nc*info_->nvirt()));
  }

  data_ = fill_block<2,double>(f, {0,0}, blocks_);
  h1_   = fill_block<2,double>(h1, {0,0}, blocks_);

  // print out the Fock matrix and terminate when preparing for external RDM runs
  if (info_->external_rdm() == "noref") {
    stringstream ss; ss << scientific << setprecision(15);
    for (int i = 0; i != nact; ++i) {
      for (int j = 0; j != nact; ++j)
        ss << f->element(j+nclosed+ncore, i+nclosed+ncore) << " ";
      ss << endl;
    }
    ofstream fs("FOCKMAT");
    fs << ss.str();
    throw Termination("Fock matrix has been dumped in FOCKMAT");
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// explict instantiation at the end of the file
template class bagel::SMITH::K2ext<double>;
template class bagel::SMITH::K2ext<complex<double>>;
template class bagel::SMITH::MOFock<double>;
template class bagel::SMITH::MOFock<complex<double>>;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

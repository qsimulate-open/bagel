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
  data_ = make_shared<TATensor<DataType,4>>(blocks_, false);
  init();
}


namespace TiledArray { namespace detail {
class DFPmap : public Pmap {
  protected:
    vector<int> ainfo_;
  public:
    DFPmap(World& world, const std::vector<int>& ainfo, const Pmap::size_type dim23) : Pmap(world, ainfo.size()*dim23), ainfo_(ainfo) {
      Pmap::local_.reserve(std::count(ainfo.begin(), ainfo.end(), Pmap::rank_) * dim23);
      for (int j = 0, cnt = 0; j != dim23; ++j)
        for (int i = 0; i != ainfo.size(); ++i, ++cnt)
          if (ainfo[i] == Pmap::rank_)
            Pmap::local_.push_back(cnt);
    }
    virtual ~DFPmap() { }
    virtual Pmap::size_type owner(const Pmap::size_type tile) const { return ainfo_[tile % ainfo_.size()]; }
    virtual bool is_local(const Pmap::size_type tile) const { return DFPmap::owner(tile) == Pmap::rank_; }
};
}}


template<>
void K2ext<complex<double>>::init() {

  vector<int> ainfo;
  IndexRange aux;
  shared_ptr<const DFDist> df = info_->geom()->df();
  const vector<pair<size_t, size_t>> atable = df->adist_now()->atable();
  for (int n = 0; n != atable.size(); ++n) {
    IndexRange a("o", atable[n].second, info_->maxtile(), aux.nblock(), aux.size());
    for (int j = 0; j != a.nblock(); ++j)
      ainfo.push_back(n);
    aux.merge(a);
  }
  auto pmap = make_shared<TiledArray::detail::DFPmap>(madness::World::get_default(), ainfo, blocks_[0].nblock()*blocks_[1].nblock());

  // Coulomb
  data_->fill_local(0.0);
  {
#ifdef HAVE_MKL_H
    mkl_set_num_threads(info_->num_threads());
#endif
    TATensor<complex<double>,3> ext(vector<IndexRange>{aux, blocks_[0], blocks_[1]}, false, pmap);
    for (auto& i0 : blocks_[0]) {
      shared_ptr<const ZMatrix> i0coeff = coeff_->slice_copy(i0.offset(), i0.offset()+i0.size());
      list<shared_ptr<RelDFHalf>> half;
      tie(half, ignore) = RelMOFile::compute_half(info_->geom(), i0coeff, /*gaunt*/false, /*breit*/false);

      for (auto& i1 : blocks_[1]) {
        shared_ptr<const ZMatrix> i1coeff = coeff_->slice_copy(i1.offset(), i1.offset()+i1.size());
        shared_ptr<const ListRelDFFull> dffull = RelMOFile::compute_full(i1coeff, half, true);
        for (auto& a : aux) {
          auto it = ext.get_local(vector<Index>{a, i0, i1});
          if (it.first)
            ext.init_tile(it.second, dffull->data().front()->get_block(a.offset(), a.size(), 0, i0.size(), 0, i1.size()));
        }
      }
    }
#ifdef HAVE_MKL_H
    mkl_set_num_threads(1);
#endif
    TATensor<complex<double>,3> ext2(vector<IndexRange>{aux, blocks_[0], blocks_[1]}, true);
    ext2("o4,o0,o1") += ext("o4,o0,o1").conj();
    (*data_)("o0,o1,o2,o3") += ext2("o4,o0,o1") * ext2("o4,o2,o3");
  }

  if (info_->gaunt()) {
    using MapType = map<int, shared_ptr<TATensor<complex<double>,3>>>;
    MapType ext, ext2;

#ifdef HAVE_MKL_H
    mkl_set_num_threads(info_->num_threads());
#endif
    {
      auto temp = make_shared<TATensor<complex<double>,3>>(vector<IndexRange>{aux, blocks_[0], blocks_[1]}, false, pmap);
      ext = MapType{{Comp::X, temp}, {Comp::Y, temp->clone(pmap)}, {Comp::Z, temp->clone(pmap)}};

      if (info_->breit())
        ext2 = MapType{{Comp::X, temp->clone(pmap)}, {Comp::Y, temp->clone(pmap)}, {Comp::Z, temp->clone(pmap)}};

      for (auto& i0 : blocks_[0]) {
        shared_ptr<const ZMatrix> i0coeff = coeff_->slice_copy(i0.offset(), i0.offset()+i0.size());
        list<shared_ptr<RelDFHalf>> half, half2;
        tie(half, half2) = RelMOFile::compute_half(info_->geom(), i0coeff, true, info_->breit());

        for (auto& i1 : blocks_[1]) {
          shared_ptr<const ZMatrix> i1coeff = coeff_->slice_copy(i1.offset(), i1.offset()+i1.size());
          auto init_local = [&](list<shared_ptr<RelDFFull>> li, MapType ee) {
            for (auto& c : li)
              for (auto& a : aux) {
                auto it = ee.at(c->alpha_comp())->get_local(vector<Index>{a, i0, i1});
                if (it.first)
                  ee.at(c->alpha_comp())->init_tile(it.second, c->get_block(a.offset(), a.size(), 0, i0.size(), 0, i1.size()));
              }
          };
          init_local(RelMOFile::compute_full(i1coeff, half, true)->data(), ext);
          if (info_->breit())
            init_local(RelMOFile::compute_full(i1coeff, half2, false)->data(), ext2);
        }
      }
      if (!info_->breit())
        ext2 = ext;
    }
#ifdef HAVE_MKL_H
    mkl_set_num_threads(1);
#endif
    auto make_conj = [](shared_ptr<const TATensor<complex<double>,3>> o) {
      auto out = o->clone();
      out->fill_local(0.0);
      (*out)("o4,o0,o1") += (*o)("o4,o0,o1").conj();
      return out;
    };
    for (auto& i : ext)
      i.second = make_conj(i.second);
    if (info_->breit())
      for (auto& i : ext2)
        i.second = make_conj(i.second);

    const double scale = info_->breit() ? -0.25 /*we explicitly symmetrize*/ : -1.0;
    (*data_)("o0,o1,o2,o3") += (*ext[Comp::X])("o4,o0,o1") * ((*ext2[Comp::X])("o4,o2,o3") * scale);
    (*data_)("o0,o1,o2,o3") += (*ext[Comp::Y])("o4,o0,o1") * ((*ext2[Comp::Y])("o4,o2,o3") * scale);
    (*data_)("o0,o1,o2,o3") += (*ext[Comp::Z])("o4,o0,o1") * ((*ext2[Comp::Z])("o4,o2,o3") * scale);
    if (info_->breit()) {
      (*data_)("o0,o1,o2,o3") += (*ext2[Comp::X])("o4,o0,o1") * ((*ext[Comp::X])("o4,o2,o3") * scale);
      (*data_)("o0,o1,o2,o3") += (*ext2[Comp::Y])("o4,o0,o1") * ((*ext[Comp::Y])("o4,o2,o3") * scale);
      (*data_)("o0,o1,o2,o3") += (*ext2[Comp::Z])("o4,o0,o1") * ((*ext[Comp::Z])("o4,o2,o3") * scale);
    }
  }

}



template<>
void K2ext<double>::init() {
  shared_ptr<const DFDist> df = info_->geom()->df();

  // It is the easiest to do integral transformation for each blocks.
  assert(blocks_.size() == 4);
  map<size_t, shared_ptr<DFFullDist>> dflist;
  // AO dimension
  assert(df->nbasis0() == df->nbasis1());

  vector<int> ainfo;
  IndexRange aux;
  const vector<pair<size_t, size_t>> atable = df->adist_now()->atable();
  for (int n = 0; n != atable.size(); ++n) {
    IndexRange a("o", atable[n].second, info_->maxtile(), aux.nblock(), aux.size());
    for (int j = 0; j != a.nblock(); ++j)
      ainfo.push_back(n);
    aux.merge(a);
  }
  auto pmap = make_shared<TiledArray::detail::DFPmap>(madness::World::get_default(), ainfo, blocks_[0].nblock()*blocks_[1].nblock());
  TATensor<double,3> ext(vector<IndexRange>{aux, blocks_[0], blocks_[1]}, false, pmap);

#ifdef HAVE_MKL_H
  mkl_set_num_threads(info_->num_threads());
#endif
  // occ loop
  for (auto& i0 : blocks_[0]) {
    shared_ptr<DFHalfDist> df_half = df->compute_half_transform(coeff_->slice(i0.offset(), i0.offset()+i0.size()))->apply_J();
    // virtual loop
    for (auto& i1 : blocks_[1]) {
      shared_ptr<DFFullDist> df_full = df_half->compute_second_transform(coeff_->slice(i1.offset(), i1.offset()+i1.size()));
      for (auto& a : aux) {
        auto it = ext.get_local(vector<Index>{a, i0, i1});
        if (it.first)
          ext.init_tile(it.second, df_full->get_block(a.offset(), a.size(), 0, i0.size(), 0, i1.size()));
      }
    }
  }
  data_->fill_local(0.0);
#ifdef HAVE_MKL_H
  mkl_set_num_threads(1);
#endif
  TATensor<double,3> ext2(vector<IndexRange>{aux, blocks_[0], blocks_[1]}, true);
  ext2("o4,o0,o1") += ext("o4,o0,o1");
  (*data_)("o0,o1,o2,o3") += ext2("o4,o0,o1") * ext2("o4,o2,o3");

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
    eig_ = VectorB(n);
    for (int i = 0; i != n; ++i) {
      assert(fabs(imag(eig[i])) < 1.0e-10);
      eig_(i) = real(eig[i]);
    }
    // move ncore to first
    if (ncore) {
      VectorB tmp(n);
      copy_n(eig_.data(), ncore, tmp.data());
      copy_n(eig_.data()+ncore+nclosed, ncore, tmp.data()+ncore);
      copy_n(eig_.data()+ncore, nclosed, tmp.data()+ncore*2);
      copy_n(eig_.data()+ncore*2+nclosed, nclosed, tmp.data()+ncore*2+nclosed);
      copy_n(eig_.data()+ncore*2+nclosed*2, nact*2+nvirt*2, tmp.data()+ncore*2+nclosed*2);
      eig_ = tmp;
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
    eig_ = VectorB(n);
    copy_n(eig.get(), n, eig_.data());
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

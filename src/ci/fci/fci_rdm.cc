//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_rdm.cc
// Copyright (C) 2011 Toru Shiozaki
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


#include <src/ci/fci/fci.h>
#include <src/util/prim_op.h>
#include <src/util/math/algo.h>
#include <src/wfn/rdm.h>

using namespace std;
using namespace bagel;

FCI_bare::FCI_bare(shared_ptr<const CIWfn> ci) {
  print_thresh_ = 1.0e-8;
  nelea_ = ci->det()->nelea();
  neleb_ = ci->det()->neleb();
  ncore_ = ci->ncore();
  norb_  = ci->nact();
  nstate_ = ci->nstates();
  energy_ = ci->energies();
  cc_ = ci->civectors()->copy();
  det_ = ci->det();
  rdm1_ = make_shared<VecRDM<1>>();
  rdm2_ = make_shared<VecRDM<2>>();
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  FCI::compute_rdm12_from_civec(shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {

  // since we consider here number conserving operators...
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<Dvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = make_shared<Dvec>(cket->det(), norb_*norb_);
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }

  return compute_rdm12_last_step(dbra, dket, cbra);
}


// computes 3 and 4RDM
tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> FCI::rdm34(const int ist, const int jst) const {
  auto rdm3 = make_shared<RDM<3>>(norb_);
  auto rdm4 = make_shared<RDM<4>>(norb_);

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);

  shared_ptr<Civec> cbra = cc_->data(ist);
  shared_ptr<Civec> cket = cc_->data(jst);

  // first make <I|E_ij|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<Dvec> dket = dbra;
  if (cbra != cket) {
    dket = dbra->clone();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  }

  // second make <J|E_kl|I><I|E_ij|0> - delta_li <J|E_kj|0>
  auto make_evec = [this](shared_ptr<Dvec> d, shared_ptr<Matrix> e, shared_ptr<Dvec> tmp) {
    int ijkl = 0;
    int ij = 0;
    for (auto iter = d->dvec().begin(); iter != d->dvec().end(); ++iter, ++ij) {
      if (ij % mpi__->size() == mpi__->rank()) {
        const int j = ij/norb_;
        const int i = ij-j*norb_;
        tmp->zero();
        sigma_2a1(*iter, tmp);
        sigma_2a2(*iter, tmp);
        int kl = 0;
        for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijkl, ++kl) {
          copy_n((*t)->data(), e->ndim(), e->element_ptr(0,ijkl));
          const int l = kl/norb_;
          const int k = kl-l*norb_;
          if (l == i)
            blas::ax_plus_y_n(-1.0, d->data(k+j*norb_)->data(), e->ndim(), e->element_ptr(0,ijkl));
        }
      } else {
        ijkl += tmp->dvec().size();
      }
    }
    e->allreduce();
  };
  auto ebra = make_shared<Matrix>(cbra->det()->size(), norb_*norb_*norb_*norb_);
  auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
  make_evec(dbra, ebra, tmp);

  shared_ptr<Matrix> eket = ebra;
  if (cbra != cket) {
    eket = ebra->clone();
    make_evec(dket, eket, tmp);
  }

  // size of the RI space
  auto dbram = make_shared<Matrix>(dbra->det()->size(), norb_*norb_);
  copy_n(dbra->data(0)->data(), dbram->size(), dbram->data());

  // first form <0|E_mn|I><I|E_ij,kl|0>
  {
    auto tmp3 = make_shared<Matrix>(*dbram % *eket);
    sort_indices<1,0,2,0,1,1,1>(tmp3->data(), rdm3->data(), norb_, norb_, norb_*norb_*norb_*norb_);

    // then perform Eq. 49 of JCP 89 5803 (Werner's MRCI paper)
    // we assume that rdm2_[ist] is set
    for (int i0 = 0; i0 != norb_; ++i0)
      for (int i1 = 0; i1 != norb_; ++i1)
        for (int i2 = 0; i2 != norb_; ++i2)
          for (int i3 = 0; i3 != norb_; ++i3) {
            blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0, i2, i1, i0), norb_, rdm3->element_ptr(0, i3, i3, i2, i1, i0));
            blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0, i0, i3, i2), norb_, rdm3->element_ptr(0, i1, i3, i2, i1, i0));
          }
  }

  // 4RDM <0|E_ij,kl|I><I|E_mn,op|0>
  {
    auto tmp4 = make_shared<Matrix>(*ebra % *eket);
    sort_indices<1,0,3,2,4,0,1,1,1>(tmp4->data(), rdm4->data(), norb_, norb_, norb_, norb_, norb_*norb_*norb_*norb_);
    for (int l = 0; l != norb_; ++l)
      for (int k = 0; k != norb_; ++k)
        for (int j = 0; j != norb_; ++j)
          for (int b = 0; b != norb_; ++b) {
            blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,0,0,k,b,l), norb_*norb_*norb_, rdm4->element_ptr(0,0,0,j,j,k,b,l));
            blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,0,0,l,b,k), norb_*norb_*norb_, rdm4->element_ptr(0,0,0,j,b,k,j,l));
            for (int i = 0; i != norb_; ++i) {
              blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0,k,b,l), norb_, rdm4->element_ptr(0,i,b,j,i,k,j,l));
              blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0,l,b,k), norb_, rdm4->element_ptr(0,i,b,j,j,k,i,l));
              for (int d = 0; d != norb_; ++d) {
                blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,k,b,j,d,l), norb_, rdm4->element_ptr(0,i,b,j,i,k,d,l));
                blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,l,b,j,d,k), norb_, rdm4->element_ptr(0,i,b,j,d,k,i,l));
              }
            }
          }
  }

  cc_->set_det(det_);

  return make_tuple(rdm3, rdm4);
}


#if 0
// computes 3 and 4RDM
// TODO duplicate code to be cleaned up
tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<3>>> FCI::rdm34f(const int ist, const int jst, shared_ptr<const Matrix> fock) const {
  auto rdm3 = make_shared<RDM<3>>(norb_);
  auto frdm4 = make_shared<RDM<3>>(norb_);

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);

  shared_ptr<Civec> cbra = cc_->data(ist);
  shared_ptr<Civec> cket = cc_->data(jst);

  // first make <I|E_ij|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<Dvec> dket = dbra;
  if (cbra != cket) {
    dket = dbra->clone();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  }

  // second make <J|E_kl|I><I|E_ij|0> - delta_li <J|E_kj|0>
  auto make_evec = [this](shared_ptr<Dvec> d, shared_ptr<Dvec> e, shared_ptr<Dvec> tmp) {
    int ijkl = 0;
    int ij = 0;
    for (auto iter = d->dvec().begin(); iter != d->dvec().end(); ++iter, ++ij) {
      const int j = ij/norb_;
      const int i = ij-j*norb_;
      tmp->zero();
      sigma_2a1(*iter, tmp);
      sigma_2a2(*iter, tmp);
      int kl = 0;
      for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijkl, ++kl) {
        *e->data(ijkl) = **t;
        const int l = kl/norb_;
        const int k = kl-l*norb_;
        if (l == i) *e->data(ijkl) -= *d->data(k+j*norb_);
      }
    }
  };
  auto ebra = make_shared<Dvec>(cbra->det(), norb_*norb_*norb_*norb_);
  auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
  make_evec(dbra, ebra, tmp);

  shared_ptr<Dvec> eket = ebra;
  if (cbra != cket) {
    eket = ebra->clone();
    make_evec(dket, eket, tmp);
  }

  // size of the RI space
  const size_t nri = ebra->lena() * ebra->lenb();
  assert(nri == dbra->lena()*dbra->lenb());

  // first form <0|E_mn|I><I|E_ij,kl|0>
  {
    auto tmp3 = make_shared<RDM<3>>(norb_);
    auto tmp3v = group(group(*tmp3,2,6),0,2);
    contract(1.0, group(*dbra,0,2), {0,1}, group(*eket,0,2), {0,2}, 0.0, tmp3v, {1,2});
    sort_indices<1,0,2,0,1,1,1>(tmp3->data(), rdm3->data(), norb_, norb_, norb_*norb_*norb_*norb_);

    // then perform Eq. 49 of JCP 89 5803 (Werner's MRCI paper)
    // we assume that rdm2_[ist] is set
    for (int i0 = 0; i0 != norb_; ++i0)
      for (int i1 = 0; i1 != norb_; ++i1)
        for (int i2 = 0; i2 != norb_; ++i2)
          for (int i3 = 0; i3 != norb_; ++i3) {
            blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0, i2, i1, i0), norb_, rdm3->element_ptr(0, i3, i3, i2, i1, i0));
            blas::ax_plus_y_n(-1.0, rdm2_->at(ist, jst)->element_ptr(0, i0, i3, i2), norb_, rdm3->element_ptr(0, i1, i3, i2, i1, i0));
          }
  }

  // 4RDM <0|E_ij,kl|I><I|E_mn,op|0>
  // now make target:  [0|E_ik,jl,op|0]  =  <0|E_ik,jl|K>[K|E_op|0] - f_nk(<0|E_in,jl,op|0>+d_ol<0|E_in,jp|0>)
  //                                                                - f_nl(<0|E_ik,jn,op|0>+d_ok<0|E_ip,jn|0>)
  //                                                                - d_ok [0|E_ip,jl|0]  - d_ol [0|E_ik,jp|0]
  {
    auto rdm2 = make_shared<RDM<2>>(norb_);
    dgemv_("T", nri, eket->ij(), 1.0, eket->data(), nri, cbra->data(), 1, 0.0, rdm2->data(), 1);

    // [0|E_ip,jl|0]
    auto frdm3 = make_shared<RDM<2>>(norb_);
    auto rdm3view = group(group(*rdm3, 4,6), 0,4);
    auto frdm3view = group(*frdm3, 0,4);
    contract(1.0, rdm3view, {0,1}, group(*fock, 0,2), {1}, 0.0, frdm3view, {0});

    // <0|E_ip,jn|0>f_nl (in this order)
    auto prdm2 = make_shared<RDM<2>>(norb_);
    auto prdm2v = group(*prdm2, 0,3);
    contract(1.0, group(*rdm2, 0,3), {0,1}, *fock, {1,2}, 0.0, prdm2v, {0,2});

    // <0|E_ik,jl,on|0>f_np (in this order)
    auto prdm3 = make_shared<RDM<3>>(norb_);
    auto prdm3v = group(*prdm3, 0,5);
    contract(1.0, group(*rdm3, 0,5), {0,1}, *fock, {1,2}, 0.0, prdm3v, {0,2});

    // <0|E_ik,jl|K>[K|E_op|0]
    auto tmp = make_shared<RDM<3>>(norb_);
    auto tmpv = group(group(*tmp, 4,6), 0,4);
    auto feket = dket->clone();
    dgemv_("N", nri*norb_*norb_, norb_*norb_, 1.0, eket->data(), nri*norb_*norb_, fock->data(), 1, 0.0, feket->data(), 1);
    contract(1.0, group(*ebra, 0,2), {0,1}, group(*feket, 0,2), {0,2}, 0.0, tmpv, {1,2});
    sort_indices<1,0,3,2,4,0,1,1,1>(tmp->data(), frdm4->data(), norb_, norb_, norb_, norb_, norb_*norb_);

    sort_indices<2,0,1,1,1,-1,1>(prdm3->data(), frdm4->data(), norb_*norb_, norb_*norb_, norb_*norb_);
    sort_indices<0,2,1,1,1,-1,1>(prdm3->data(), frdm4->data(), norb_*norb_, norb_*norb_, norb_*norb_);

    auto prdm2t = prdm2->clone();
    sort_indices<1,0,0,1,1,1>(prdm2->data(), prdm2t->data(), norb_*norb_, norb_*norb_);
    for (int p = 0; p != norb_; ++p)
      for (int l = 0; l != norb_; ++l) {
        blas::ax_plus_y_n(-1.0, prdm2t->element_ptr(0,0,0,p), norb_*norb_*norb_, frdm4->element_ptr(0,0,0,l,l,p));
        blas::ax_plus_y_n(-1.0, frdm3->element_ptr(0,0,0,p),  norb_*norb_*norb_, frdm4->element_ptr(0,0,0,l,l,p));
        for (int k = 0; k != norb_; ++k)
          for (int j = 0; j != norb_; ++j) {
            blas::ax_plus_y_n(-1.0, prdm2->element_ptr(0,p,j,l), norb_, frdm4->element_ptr(0,k,j,l,k,p));
            blas::ax_plus_y_n(-1.0, frdm3->element_ptr(0,p,j,l), norb_, frdm4->element_ptr(0,k,j,l,k,p));
          }
      }
  }

  cc_->set_det(det_);

  return make_tuple(rdm3, frdm4);
}
#endif


// note that this does not transform internal integrals (since it is not needed in CASSCF).
pair<shared_ptr<Matrix>, VectorB> FCI::natorb_convert() {
  assert(rdm1_av_ != nullptr);
  pair<shared_ptr<Matrix>, VectorB> natorb = rdm1_av_->generate_natural_orbitals();
  update_rdms(natorb.first);
  jop_->update_1ext_ints(natorb.first);
  for (auto& i : natorb.second)
    if (i < numerical_zero__) i = 0.0;
  return natorb;
}


void FCI::update_rdms(const shared_ptr<Matrix>& coeff) {
  for (auto& i : *rdm1_)
    i.second->transform(coeff);
  for (auto& i : *rdm2_)
    i.second->transform(coeff);

  // Only when #state > 1, this is needed.
  // Actually rdm1_av_ points to the same object as rdm1_ in 1 state runs. Therefore if you do twice, you get wrong.
  if (rdm1_->size() > 1) rdm1_av_->transform(coeff);
  if (rdm2_->size() > 1) rdm2_av_->transform(coeff);
  assert(rdm1_->size() > 1 || rdm1_->at(0) == rdm1_av_);
  assert(rdm2_->size() > 1 || rdm2_->at(0) == rdm2_av_);
}


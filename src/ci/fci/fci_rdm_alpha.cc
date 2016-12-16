//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_rdm_alpha.cc
// Copyright (C) 2016 Toru Shiozaki
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

using namespace std;
using namespace bagel;

// This file is devoted for computing the so-called "alpha" density matrices
// <iaja>, <kl,iaja>, <kl,mn,iaja>, <kl,mn,op,iaja>
// the formulas are more or less the same.

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> FCI::rdm12_alpha(const int ist, const int jst) const {
  if (det_->compress()) {
    auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
    cc_->set_det(detex);
  }

  shared_ptr<Civec> cbra = cc_->data(ist);
  shared_ptr<Civec> cket = cc_->data(jst);

  // ij applied to bra
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  // kala applied to ket
  auto dket = make_shared<Dvec>(cket->det(), norb_*norb_);
  sigma_2a1(cket, dket);

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  tie(rdm1,rdm2) = compute_rdm12_last_step(dbra, dket, cbra);

  cc_->set_det(det_);
  return tie(rdm1, rdm2);
}


tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> FCI::rdm34_alpha(const int ist, const int jst) const {
  shared_ptr<const RDM<1>> rdm1;
  shared_ptr<const RDM<2>> rdm2;
  tie(rdm1, rdm2) = rdm12_alpha(ist, jst);

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

  // dket is alpha
  auto dket = dbra->clone();
  sigma_2a1(cket, dket);
//sigma_2a2(cket, dket); ... KET side has iaja

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

  // last two indices of eket is alpha
  auto eket = ebra->clone();
  make_evec(dket, eket, tmp);

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
    for (int i0 = 0; i0 != norb_; ++i0)
      for (int i1 = 0; i1 != norb_; ++i1)
        for (int i2 = 0; i2 != norb_; ++i2)
          for (int i3 = 0; i3 != norb_; ++i3) {
            blas::ax_plus_y_n(-1.0, rdm2->element_ptr(0, i2, i1, i0), norb_, rdm3->element_ptr(0, i3, i3, i2, i1, i0));
            for (int i4 = 0; i4 != norb_; ++i4)
              rdm3->element(i4, i1, i3, i2, i1, i0) -= rdm2->element(i3, i2, i4, i0);
          }
  }

  // 4RDM <0|E_ij,kl|I><I|E_mn,op|0>
  {
    {
      auto tmp4 = make_shared<RDM<4>>(norb_);
      dgemm_("T", "N", ebra->ij(), eket->ij(), nri, 1.0, ebra->data(), nri, eket->data(), nri, 0.0, tmp4->data(), ebra->ij());
      sort_indices<1,0,3,2,4,0,1,1,1>(tmp4->data(), rdm4->data(), norb_, norb_, norb_, norb_, norb_*norb_*norb_*norb_);
      for (int l = 0; l != norb_; ++l)
        for (int k = 0; k != norb_; ++k)
          for (int j = 0; j != norb_; ++j)
            for (int b = 0; b != norb_; ++b) {
              blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,0,0,k,b,l), norb_*norb_*norb_, rdm4->element_ptr(0,0,0,j,j,k,b,l));
              for (int i = 0; i != norb_; ++i) {
                blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,0,b,k,i,l), norb_*norb_, rdm4->element_ptr(0,0,i,j,b,k,j,l));
                blas::ax_plus_y_n(-1.0, rdm2->element_ptr(0,k,b,l), norb_, rdm4->element_ptr(0,i,b,j,i,k,j,l));
                for (int d = 0; d != norb_; ++d) {
                  rdm4->element(d,i,b,j,j,k,i,l) -= rdm2->element(b,k,d,l);
                  blas::ax_plus_y_n(-1.0, rdm3->element_ptr(0,k,b,j,d,l), norb_, rdm4->element_ptr(0,i,b,j,i,k,d,l));
                  for (int a = 0; a != norb_; ++a)
                    rdm4->element(a,i,b,j,d,k,i,l) -= rdm3->element(b,j,d,k,a,l);
                }
              }
            }
    }
  }

  cc_->set_det(det_);
  return make_tuple(rdm3, rdm4);
}

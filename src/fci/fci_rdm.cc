//
// BAGEL - Parallel electron correlation program.
// Filename: fci_rdm.cc
// Copyright (C) 2011 Toru Shiozaki
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


#include <src/fci/fci.h>
#include <src/smith/prim_op.h>
#include <src/math/algo.h>
#include <src/wfn/rdm.h>

using namespace std;
using namespace bagel;

void FCI::compute_rdm12() {
  // Needs initialization here because we use daxpy.
  // For nstate_ == 1, rdm1_av_ = rdm1_[0].
  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(norb_);
    rdm2_av_ = make_shared<RDM<2>>(norb_);
  }
  if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }
  // we need expanded lists
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);

  for (int i = 0; i != nstate_; ++i) compute_rdm12(i);

  cc_->set_det(det_);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  FCI::compute_rdm12_last_step(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Civec> cibra) const {

  const int nri = dbra->lena()*dbra->lenb();
  const int ij  = norb_*norb_;

  if (nri != dket->lena()*dket->lenb())
    throw logic_error("FCI::compute_rdm12_last_step called with inconsistent RI spaces");

  // 1RDM
  // c^dagger <I|\hat{E}|0>
  auto rdm1 = make_shared<RDM<1>>(norb_);
  dgemv_("T", nri, ij, 1.0, dket->data(0)->data(), nri, cibra->data(), 1, 0.0, rdm1->data(), 1);
  // 2RDM
  // \sum_I <0|\hat{E}|I> <I|\hat{E}|0>
  auto rdm2 = make_shared<RDM<2>>(norb_);
  dgemm_("T", "N", ij, ij, nri, 1.0, dbra->data(0)->data(), nri, dket->data(0)->data(), nri, 0.0, rdm2->data(), ij);

  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb_*norb_]);
  for (int i = 0; i != norb_; ++i) {
    for (int k = 0; k != norb_; ++k) {
      copy_n(&rdm2->element(0,0,k,i), norb_*norb_, buf.get());
      blas::transpose(buf.get(), norb_, norb_, &rdm2->element(0,0,k,i)); // sorting with stride 1 as norb_ is small
    }
  }

  // put in diagonal into 2RDM
  // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
  for (int i = 0; i != norb_; ++i)
    for (int k = 0; k != norb_; ++k)
      for (int j = 0; j != norb_; ++j)
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  return tie(rdm1, rdm2);
}


shared_ptr<Dvec> FCI::rdm1deriv() const {

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<Civec> cbra = cc_->data(0);

  // 1RDM ci derivative
  // <I|E_ij|0>

  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  return dbra;
}


shared_ptr<Dvec> FCI::rdm2deriv() const {

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<Civec> cbra = cc_->data(0);

  // make  <I|E_ij|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  // second make <J|E_kl|I><I|E_ij|0> - delta_li <J|E_kj|0>
  auto ebra = make_shared<Dvec>(cbra->det(), norb_*norb_*norb_*norb_);
  auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
  int ijkl = 0;
  int ij = 0;
  for (auto iter = dbra->dvec().begin(); iter != dbra->dvec().end(); ++iter, ++ij) {
    const int j = ij/norb_;
    const int i = ij-j*norb_;
    tmp->zero();
    sigma_2a1(*iter, tmp);
    sigma_2a2(*iter, tmp);
    int kl = 0;
    for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijkl, ++kl) {
      *ebra->data(ijkl) = **t;
      const int l = kl/norb_;
      const int k = kl-l*norb_;
      if (l == i) *ebra->data(ijkl) -= *dbra->data(k+norb_*j);
    }
  }
  return ebra;
}


shared_ptr<Dvec> FCI::rdm3deriv() const {
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<Civec> cbra = cc_->data(0);

  // first make <I|i+j|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  // second make <J|k+i+jl|0> = <J|k+l|I><I|i+j|0> - delta_li <J|k+j|0>
  auto ebra = make_shared<Dvec>(cbra->det(), norb_*norb_*norb_*norb_);
  {
    auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
    int ijkl = 0;
    int ij = 0;
    for (auto iter = dbra->dvec().begin(); iter != dbra->dvec().end(); ++iter, ++ij) {
      const int j = ij/norb_;
      const int i = ij-j*norb_;
      tmp->zero();
      sigma_2a1(*iter, tmp);
      sigma_2a2(*iter, tmp);
      int kl = 0;
      for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijkl, ++kl) {
        *ebra->data(ijkl) = **t;
        const int l = kl/norb_;
        const int k = kl-l*norb_;
        if (l == i) *ebra->data(ijkl) -= *dbra->data(k+norb_*j);
      }
    }
  }

  // now make target  <K|m+k+i+jln|0>  =  <K|m+n|J><J|k+i+jl|0> - delta_nk<K|m+i+jl|0> - delta_ni<K|k+m+jl|0>
  auto fbra = make_shared<Dvec>(cbra->det(), norb_*norb_*norb_*norb_*norb_*norb_);
  {
    auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
    int ijklmn = 0;
    int ijkl = 0;
    for (auto iter = ebra->dvec().begin(); iter != ebra->dvec().end(); ++iter, ++ijkl) {
      const int j = ijkl/(norb_*norb_*norb_);
      const int i = ijkl/(norb_*norb_)-j*norb_;
      const int l = ijkl/norb_-i*norb_-j*norb_*norb_;
      const int k = ijkl-l*norb_-i*norb_*norb_-j*norb_*norb_*norb_;
      tmp->zero();
      sigma_2a1(*iter, tmp);
      sigma_2a2(*iter, tmp);
      int mn = 0;
      for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijklmn, ++mn) {
        const int n = mn/norb_;
        const int m = mn-n*norb_;
        *fbra->data(ijklmn) = **t;
        if (n == k) *fbra->data(ijklmn) -= *ebra->data(m+norb_*(l+norb_*(i+norb_*(j))));
        if (n == i) *fbra->data(ijklmn) -= *ebra->data(k+norb_*(l+norb_*(m+norb_*(j))));
      }
    }
  }

  return fbra;
}


shared_ptr<Dvec> FCI::rdm4deriv() const {
  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);
  shared_ptr<Civec> cbra = cc_->data(0);

  // first make <I|i+j|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  // second make <J|k+i+jl|0> = <J|k+l|I><I|i+j|0> - delta_li <J|k+j|0>
  auto ebra = make_shared<Dvec>(cbra->det(), norb_*norb_*norb_*norb_);
  {
    auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
    int ijkl = 0;
    int ij = 0;
    for (auto iter = dbra->dvec().begin(); iter != dbra->dvec().end(); ++iter, ++ij) {
      const int j = ij/norb_;
      const int i = ij-j*norb_;
      tmp->zero();
      sigma_2a1(*iter, tmp);
      sigma_2a2(*iter, tmp);
      int kl = 0;
      for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijkl, ++kl) {
        *ebra->data(ijkl) = **t;
        const int l = kl/norb_;
        const int k = kl-l*norb_;
        if (l == i) *ebra->data(ijkl) -= *dbra->data(k+norb_*j);
      }
    }
  }

  // third make <K|m+k+i+jln|0>  =  <K|m+n|J><J|k+i+jl|0> - delta_nk<K|m+i+jl|0> - delta_ni<K|k+m+jl|0>
  auto fbra = make_shared<Dvec>(cbra->det(), norb_*norb_*norb_*norb_*norb_*norb_);
  {
    auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
    int ijklmn = 0;
    int ijkl = 0;
    for (auto iter = ebra->dvec().begin(); iter != ebra->dvec().end(); ++iter, ++ijkl) {
      const int j = ijkl/(norb_*norb_*norb_);
      const int i = ijkl/(norb_*norb_)-j*norb_;
      const int l = ijkl/norb_-i*norb_-j*norb_*norb_;
      const int k = ijkl-l*norb_-i*norb_*norb_-j*norb_*norb_*norb_;
      tmp->zero();
      sigma_2a1(*iter, tmp);
      sigma_2a2(*iter, tmp);
      int mn = 0;
      for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijklmn, ++mn) {
        const int n = mn/norb_;
        const int m = mn-n*norb_;
        *fbra->data(ijklmn) = **t;
        if (n == k) *fbra->data(ijklmn) -= *ebra->data(m+norb_*(l+norb_*(i+norb_*(j))));
        if (n == i) *fbra->data(ijklmn) -= *ebra->data(k+norb_*(l+norb_*(m+norb_*(j))));
      }
    }
  }

  // now make target:  <L|o+m+k+i+jlnp|0>  =  <L|o+p|K><K|m+k+i+jln|0> - delta_pm<L|o+k+i+jln|0> - delta_pk<L|m+o+i+jln|0> - delta_pi<L|m+k+o+jln|0>
  auto gbra = make_shared<Dvec>(cbra->det(), norb_*norb_*norb_*norb_*norb_*norb_*norb_*norb_);
  {
    auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
    int ijklmnop = 0;
    int ijklmn = 0;
    for (auto iter = fbra->dvec().begin(); iter != fbra->dvec().end(); ++iter, ++ijklmn) {
      const int j = ijklmn/(norb_*norb_*norb_*norb_*norb_);
      const int i = ijklmn/(norb_*norb_*norb_*norb_)-j*norb_;
      const int l = ijklmn/(norb_*norb_*norb_)-i*norb_-j*norb_*norb_;
      const int k = ijklmn/(norb_*norb_)-l*norb_-i*norb_*norb_-j*norb_*norb_*norb_;
      const int n = ijklmn/norb_-k*norb_-l*norb_*norb_-i*norb_*norb_*norb_-j*norb_*norb_*norb_*norb_;
      const int m = ijklmn-n*norb_-k*norb_*norb_-l*norb_*norb_*norb_-i*norb_*norb_*norb_*norb_-j*norb_*norb_*norb_*norb_*norb_;
      tmp->zero();
      sigma_2a1(*iter, tmp);
      sigma_2a2(*iter, tmp);
      int op = 0;
      for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijklmnop, ++op) {
        const int p = op/norb_;
        const int o = op-p*norb_;
        *gbra->data(ijklmnop) = **t;
        if (p == m) *gbra->data(ijklmnop) -= *fbra->data(o+norb_*(n+norb_*(k+norb_*(l+norb_*(i+norb_*(j))))));
        if (p == k) *gbra->data(ijklmnop) -= *fbra->data(m+norb_*(n+norb_*(o+norb_*(l+norb_*(i+norb_*(j))))));
        if (p == i) *gbra->data(ijklmnop) -= *fbra->data(m+norb_*(n+norb_*(k+norb_*(l+norb_*(o+norb_*(j))))));
      }
    }
  }

  return gbra;
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  FCI::compute_rdm12_from_civec(shared_ptr<const Civec> cbra, shared_ptr<const Civec> cket) const {

  // since we consider here number conserving operators...
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  shared_ptr<Dvec> dket;
  // if bra and ket vectors are different, we need to form Sigma for ket as well.
  if (cbra != cket) {
    dket = make_shared<Dvec>(cket->det(), norb_*norb_);
    dket->zero();
    sigma_2a1(cket, dket);
    sigma_2a2(cket, dket);
  } else {
    dket = dbra;
  }

  return compute_rdm12_last_step(dbra, dket, cbra);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  FCI::compute_rdm12_av_from_dvec(shared_ptr<const Dvec> dbra, shared_ptr<const Dvec> dket, shared_ptr<const Determinants> o) const {

  if (o != nullptr) {
    dbra->set_det(o);
    dket->set_det(o);
  }

  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);
  rdm1->zero();
  rdm2->zero();

  assert(dbra->ij() == dket->ij() && dbra->det() == dket->det());

  for (int i = 0; i != dbra->ij(); ++i) {
    shared_ptr<RDM<1>> r1;
    shared_ptr<RDM<2>> r2;
    tie(r1, r2) = compute_rdm12_from_civec(dbra->data(i), dket->data(i));
    rdm1->ax_plus_y(weight_[i], r1);
    rdm2->ax_plus_y(weight_[i], r2);
  }

  if (o != nullptr) {
    dbra->set_det(det_);
    dket->set_det(det_);
  }

  return tie(rdm1, rdm2);
}


void FCI::compute_rdm12(const int ist) {
  shared_ptr<Civec> cc = cc_->data(ist);

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  tie(rdm1, rdm2) = compute_rdm12_from_civec(cc, cc);

  // setting to private members.
  rdm1_[ist] = rdm1;
  rdm2_[ist] = rdm2;
  if (nstate_ != 1) {
    rdm1_av_->ax_plus_y(weight_[ist], rdm1);
    rdm2_av_->ax_plus_y(weight_[ist], rdm2);
  } else {
    rdm1_av_ = rdm1;
    rdm2_av_ = rdm2;
  }

}

// computes 3 and 4RDM
tuple<shared_ptr<RDM<3>>, shared_ptr<RDM<4>>> FCI::compute_rdm34(const int ist) const {
  auto rdm3 = make_shared<RDM<3>>(norb_);
  auto rdm4 = make_shared<RDM<4>>(norb_);

  auto detex = make_shared<Determinants>(norb_, nelea_, neleb_, false, /*mute=*/true);
  cc_->set_det(detex);

  shared_ptr<Civec> cbra = cc_->data(ist);

  // first make <I|E_ij|0>
  auto dbra = make_shared<Dvec>(cbra->det(), norb_*norb_);
  dbra->zero();
  sigma_2a1(cbra, dbra);
  sigma_2a2(cbra, dbra);

  // second make <J|E_kl|I><I|E_ij|0> - delta_li <J|E_kj|0>
  auto ebra = make_shared<Dvec>(cbra->det(), norb_*norb_*norb_*norb_);
  auto tmp = make_shared<Dvec>(cbra->det(), norb_*norb_);
  int ijkl = 0;
  int ij = 0;
  for (auto iter = dbra->dvec().begin(); iter != dbra->dvec().end(); ++iter, ++ij) {
    const int j = ij/norb_;
    const int i = ij-j*norb_;
    tmp->zero();
    sigma_2a1(*iter, tmp);
    sigma_2a2(*iter, tmp);
    int kl = 0;
    for (auto t = tmp->dvec().begin(); t != tmp->dvec().end(); ++t, ++ijkl, ++kl) {
      *ebra->data(ijkl) = **t;
      const int l = kl/norb_;
      const int k = kl-l*norb_;
      if (l == i) *ebra->data(ijkl) -= *dbra->data(k+j*norb_);
    }
  }

  // size of the RI space
  const size_t nri = ebra->lena() * ebra->lenb();
  assert(nri == dbra->lena()*dbra->lenb());

  // first form <0|E_ij,kl|I><I|E_mn|0>
  {
    auto tmp3 = make_shared<RDM<3>>(norb_);
    dgemm_("T", "N", dbra->ij(), ebra->ij(), nri, 1.0, dbra->data(), nri, ebra->data(), nri, 0.0, tmp3->data(), dbra->ij());

    // then perform Eq. 49 of JCP 89 5803 (Werner's MRCI paper)
    // we assume that rdm2_[ist] is set
    for (int i0 = 0; i0 != norb_; ++i0) {
      for (int i1 = 0; i1 != norb_; ++i1) {
        for (int i2 = 0; i2 != norb_; ++i2) {
          for (int i3 = 0; i3 != norb_; ++i3) {
            // i4 and i5 correspond to m and n (they should be transposed here)
            for (int i5 = 0; i5 != norb_; ++i5) {
              for (int i4 = 0; i4 != norb_; ++i4) {
                rdm3->element(i5, i4, i3, i2, i1, i0) = tmp3->element(i4, i5, i3, i2, i1, i0);
              }
              rdm3->element(i5, i3, i3, i2, i1, i0) -= rdm2_[ist]->element(i5, i2, i1, i0);
              rdm3->element(i5, i1, i3, i2, i1, i0) -= rdm2_[ist]->element(i3, i2, i5, i0);
            }
          }
        }
      }
    }
  }

  // 4RDM <0|E_ij,kl|I><I|E_mn,op|0>
  {
    {
      auto tmp4 = make_shared<RDM<4>>(norb_);
      dgemm_("T", "N", ebra->ij(), ebra->ij(), nri, 1.0, ebra->data(), nri, ebra->data(), nri, 0.0, tmp4->data(), ebra->ij());
      SMITH::sort_indices<1,0,3,2,4,5,6,7,0,1,1,1>(tmp4->data(), rdm4->data(), norb_, norb_, norb_, norb_, norb_, norb_, norb_, norb_);
      for (int l = 0; l != norb_; ++l)
        for (int d = 0; d != norb_; ++d)
          for (int k = 0; k != norb_; ++k)
            for (int c = 0; c != norb_; ++c)
              for (int j = 0; j != norb_; ++j)
                for (int b = 0; b != norb_; ++b)
                  for (int i = 0; i != norb_; ++i)
                    for (int a = 0; a != norb_; ++a) {
                      if (c == i && d == j) rdm4->element(a,i,b,j,c,k,d,l) -= rdm2_[ist]->element(a,k,b,l);
                      if (c == j && d == i) rdm4->element(a,i,b,j,c,k,d,l) -= rdm2_[ist]->element(a,l,b,k);
                      if (c == i)           rdm4->element(a,i,b,j,c,k,d,l) -= rdm3->element(a,k,b,j,d,l);
                      if (c == j)           rdm4->element(a,i,b,j,c,k,d,l) -= rdm3->element(a,i,b,k,d,l);
                      if (d == i)           rdm4->element(a,i,b,j,c,k,d,l) -= rdm3->element(a,l,b,j,c,k);
                      if (d == j)           rdm4->element(a,i,b,j,c,k,d,l) -= rdm3->element(a,i,b,l,c,k);
                    }
    }
  }
#if 0
  // Checking 4RDM by comparing with 3RDM
  auto debug = make_shared<RDM<3>>(*rdm3);
  cout << "printing out rdm" << endl;
  for (int l = 0; l != norb_; ++l)
    for (int d = 0; d != norb_; ++d)
      for (int k = 0; k != norb_; ++k)
        for (int c = 0; c != norb_; ++c)
          for (int j = 0; j != norb_; ++j)
            for (int b = 0; b != norb_; ++b)
    for (int i = 0; i != norb_; ++i) {
      debug->element(b,j,c,k,d,l) -= 1.0/(nelea()+neleb()-3) * rdm4->element(i,i,b,j,c,k,d,l);
//    debug->element(b,j,c,k,d,l) -= 1.0/(nelea()+neleb()-3) * rdm4->element(b,j,i,i,c,k,d,l);
//    debug->element(b,j,c,k,d,l) -= 1.0/(nelea()+neleb()-3) * rdm4->element(b,j,c,k,i,i,d,l);
//    debug->element(b,j,c,k,d,l) -= 1.0/(nelea()+neleb()-3) * rdm4->element(b,j,c,k,d,l,i,i);
    }
  debug->print(1.0e-8);
  cout << "printing out rdm - end" << endl;
#endif

  cc_->set_det(det_);

  return make_tuple(rdm3, rdm4);
}

// note that this does not transform internal integrals (since it is not needed in CASSCF).
pair<shared_ptr<Matrix>, vector<double>> FCI::natorb_convert() {
  assert(rdm1_av_ != nullptr);
  pair<shared_ptr<Matrix>, vector<double>> natorb = rdm1_av_->generate_natural_orbitals();
  update_rdms(natorb.first);
  jop_->update_1ext_ints(natorb.first);
  return natorb;
}


void FCI::update_rdms(const shared_ptr<Matrix>& coeff) {
  for (auto iter = rdm1_.begin(); iter != rdm1_.end(); ++iter)
    (*iter)->transform(coeff);
  for (auto iter = rdm2_.begin(); iter != rdm2_.end(); ++iter)
    (*iter)->transform(coeff);

  // Only when #state > 1, this is needed.
  // Actually rdm1_av_ points to the same object as rdm1_ in 1 state runs. Therefore if you do twice, you get wrong.
  if (rdm1_.size() > 1) rdm1_av_->transform(coeff);
  if (rdm2_.size() > 1) rdm2_av_->transform(coeff);
  assert(rdm1_.size() > 1 || rdm1_.front() == rdm1_av_);
  assert(rdm2_.size() > 1 || rdm2_.front() == rdm2_av_);
}

